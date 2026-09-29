"""Pressure-folder selectors: all, a value, an inclusive range, or a list."""
import math


def parse_pressure_selection(value):
    """Return None for all, or a tuple of inclusive (low, high) intervals.

    Numeric values remain supported for callers using the original float API.
    """
    if value is None or str(value).strip().lower() == 'all':
        return None
    specification = str(value).strip()
    if not specification:
        raise ValueError('Pressure selection cannot be empty')
    intervals = []
    for item in specification.split(','):
        bounds = item.strip().split(':')
        if len(bounds) not in (1, 2) or any(not x.strip() for x in bounds):
            raise ValueError(f'Invalid pressure selection {item!r}; use all, 1, 1:5, or 0.5,1,3')
        try:
            numbers = [float(x) for x in bounds]
        except ValueError as error:
            raise ValueError(f'Invalid pressure selection {item!r}; values must be numeric') from error
        if any(not math.isfinite(x) or x <= 0 for x in numbers):
            raise ValueError('Pressure values must be finite and positive (atm)')
        low, high = numbers[0], numbers[-1]
        if low > high:
            raise ValueError(f'Pressure range {item!r} is reversed; use low:high')
        intervals.append((low, high))
    return tuple(intervals)


def selected_pressure_folders(folders, selection='all'):
    """Select existing folders by numeric pressure, not lexicographic order.

    Each requested value/range must match at least one existing folder. Range
    endpoints need not exist: 1.5:3 includes existing 2 and 3 atm folders.
    """
    intervals = parse_pressure_selection(selection)
    available = []
    for folder in folders:
        if folder.name.startswith('.') or not folder.is_dir():
            continue
        try:
            pressure = float(folder.name)
        except ValueError:
            continue
        if math.isfinite(pressure) and pressure > 0:
            available.append((pressure, folder))
    available.sort(key=lambda item: (item[0], item[1].name))
    if not available:
        raise ValueError('No pressure folders found in the database')
    if intervals is None:
        return [folder for _, folder in available]
    for low, high in intervals:
        if not any(low <= pressure <= high for pressure, _ in available):
            request = str(low) if low == high else f'{low}:{high}'
            names = ', '.join(folder.name for _, folder in available)
            raise ValueError(f'No pressure folders match {request} atm. Available: {names}')
    return [folder for pressure, folder in available
            if any(low <= pressure <= high for low, high in intervals)]
