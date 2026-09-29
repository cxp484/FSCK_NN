! ONNX inference callable from Fortran. No weights or activations are reimplemented.
module nn_model
  use iso_c_binding
  use iso_fortran_env, only: real64, error_unit
  implicit none
  private
  public :: read_model, predict, close_model
  type(c_ptr), save :: handle=c_null_ptr
  interface
    integer(c_int) function onnx_open(path, model, message, n) bind(C,name='fsck_onnx_open')
      import
      character(c_char), intent(in) :: path(*)
      type(c_ptr), intent(out) :: model
      character(c_char), intent(out) :: message(*)
      integer(c_int), value :: n
    end function
    integer(c_int) function onnx_predict(model, x, k, message, n) bind(C,name='fsck_onnx_predict')
      import
      type(c_ptr), value :: model
      real(c_float), intent(in) :: x(7)
      real(c_float), intent(out) :: k(32)
      character(c_char), intent(out) :: message(*)
      integer(c_int), value :: n
    end function
    subroutine onnx_close(model) bind(C,name='fsck_onnx_close')
      import
      type(c_ptr), value :: model
    end subroutine
  end interface
contains
  subroutine check(status,message)
    integer(c_int), intent(in) :: status
    character(c_char), intent(in) :: message(:)
    integer :: i
    if (status==0) return
    write(error_unit,'(a)',advance='no') 'ONNX Runtime: '
    do i=1,size(message)
      if (message(i)==c_null_char) exit
      write(error_unit,'(a)',advance='no') message(i)
    enddo
    write(error_unit,*)
    error stop 1
  end subroutine
  subroutine read_model(filename)
    character(*), intent(in) :: filename
    character(c_char) :: message(2048)
    integer(c_int) :: status
    call close_model()
    status=onnx_open(trim(filename)//c_null_char,handle,message,int(size(message),c_int))
    call check(status,message)
  end subroutine
  subroutine predict(state,k)
    real(real64), intent(in) :: state(7)
    real(real64), intent(out) :: k(32)
    real(c_float) :: x(7), output(32)
    character(c_char) :: message(2048)
    integer(c_int) :: status
    x=real(state,c_float)
    status=onnx_predict(handle,x,output,message,int(size(message),c_int))
    call check(status,message)
    k=real(output,real64)
  end subroutine
  subroutine close_model()
    if (c_associated(handle)) call onnx_close(handle)
    handle=c_null_ptr
  end subroutine
end module
