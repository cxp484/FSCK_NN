/* C ABI used by Fortran ISO_C_BINDING. All ONNX errors become status codes. */
#include <onnxruntime_c_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const OrtApi *api;
    OrtEnv *env;
    OrtSessionOptions *options;
    OrtSession *session;
} Model;

void fsck_onnx_close(void *handle) {
    Model *m = (Model *)handle;
    if (!m) return;
    if (m->session) m->api->ReleaseSession(m->session);
    if (m->options) m->api->ReleaseSessionOptions(m->options);
    if (m->env) m->api->ReleaseEnv(m->env);
    free(m);
}
static int report(const OrtApi *api, OrtStatus *status, char *error, int size) {
    if (size > 0) snprintf(error, (size_t)size, "%s", api->GetErrorMessage(status));
    api->ReleaseStatus(status);
    return 1;
}
#define TRY(call) do { status = (call); if (status) goto fail; } while (0)

static OrtStatus *validate(Model *m, int input) {
    const OrtApi *api = m->api;
    OrtStatus *status = NULL;
    OrtTypeInfo *type = NULL;
    const OrtTensorTypeAndShapeInfo *tensor = NULL;
    OrtAllocator *allocator = NULL;
    char *name = NULL;
    size_t count = 0, rank = 0;
    int64_t shape[2];
    ONNXTensorElementDataType element;
    TRY(input ? api->SessionGetInputCount(m->session, &count) : api->SessionGetOutputCount(m->session, &count));
    if (count != 1) { status=api->CreateStatus(ORT_INVALID_ARGUMENT,"Expected one input and one output"); goto fail; }
    TRY(api->GetAllocatorWithDefaultOptions(&allocator));
    TRY(input ? api->SessionGetInputName(m->session,0,allocator,&name) : api->SessionGetOutputName(m->session,0,allocator,&name));
    if (strcmp(name,input ? "physical_inputs" : "kappa") != 0) {
        status=api->CreateStatus(ORT_INVALID_ARGUMENT,"Expected physical_inputs and kappa tensor names"); goto fail;
    }
    TRY(input ? api->SessionGetInputTypeInfo(m->session,0,&type) : api->SessionGetOutputTypeInfo(m->session,0,&type));
    TRY(api->CastTypeInfoToTensorInfo(type,&tensor));
    if (!tensor) { status=api->CreateStatus(ORT_INVALID_ARGUMENT,"Expected tensor input/output"); goto fail; }
    TRY(api->GetTensorElementType(tensor,&element));
    TRY(api->GetDimensionsCount(tensor,&rank));
    if (rank != 2 || element != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT) {
        status=api->CreateStatus(ORT_INVALID_ARGUMENT,"Expected rank-two float32 tensors"); goto fail;
    }
    TRY(api->GetDimensions(tensor,shape,2));
    if ((shape[0] != -1 && shape[0] != 1) || shape[1] != (input ? 7 : 32))
        status=api->CreateStatus(ORT_INVALID_ARGUMENT,"Expected [batch,7] input and [batch,32] output");
fail:
    if (name) allocator->Free(allocator,name);
    if (type) api->ReleaseTypeInfo(type);
    return status;
}
int fsck_onnx_open(const char *path, void **handle, char *error, int size) {
    OrtStatus *status = NULL;
    const OrtApi *api=OrtGetApiBase()->GetApi(ORT_API_VERSION);
    *handle=NULL;
    if (!api) { if (size>0) snprintf(error,(size_t)size,"Incompatible ONNX Runtime API version"); return 1; }
    Model *m=(Model *)calloc(1,sizeof(Model));
    if (!m) { if (size>0) snprintf(error,(size_t)size,"Cannot allocate ONNX model"); return 1; }
    m->api=api;
    TRY(api->CreateEnv(ORT_LOGGING_LEVEL_WARNING,"fsck",&m->env));
    TRY(api->CreateSessionOptions(&m->options));
    TRY(api->SetIntraOpNumThreads(m->options,1));
    TRY(api->CreateSession(m->env,path,m->options,&m->session));
    TRY(validate(m,1));
    TRY(validate(m,0));
    *handle=m;
    return 0;
fail:
    report(api,status,error,size);
    fsck_onnx_close(m);
    return 1;
}
int fsck_onnx_predict(void *handle, const float *x, float *k, char *error, int size) {
    Model *m=(Model *)handle;
    if (!m) { if (size>0) snprintf(error,(size_t)size,"ONNX model is not loaded"); return 1; }
    const OrtApi *api=m->api;
    OrtStatus *status=NULL;
    OrtMemoryInfo *memory=NULL;
    OrtValue *input=NULL,*output=NULL;
    OrtTensorTypeAndShapeInfo *info=NULL;
    int64_t shape[2]={1,7};
    size_t count=0;
    void *values=NULL;
    const char *inputs[]={"physical_inputs"}, *outputs[]={"kappa"};
    TRY(api->CreateCpuMemoryInfo(OrtArenaAllocator,OrtMemTypeDefault,&memory));
    TRY(api->CreateTensorWithDataAsOrtValue(memory,(void *)x,7*sizeof(float),shape,2,ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT,&input));
    const OrtValue *input_values[]={input};
    TRY(api->Run(m->session,NULL,inputs,input_values,1,outputs,1,&output));
    TRY(api->GetTensorTypeAndShape(output,&info));
    TRY(api->GetTensorShapeElementCount(info,&count));
    if (count!=32) { status=api->CreateStatus(ORT_INVALID_ARGUMENT,"Expected 32 output values"); goto fail; }
    TRY(api->GetTensorMutableData(output,&values));
    memcpy(k,values,32*sizeof(float));
fail:
    if (info) api->ReleaseTensorTypeAndShapeInfo(info);
    if (output) api->ReleaseValue(output);
    if (input) api->ReleaseValue(input);
    if (memory) api->ReleaseMemoryInfo(memory);
    if (status) return report(api,status,error,size);
    return 0;
}
