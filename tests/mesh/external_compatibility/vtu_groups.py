"""Independent VTK access to the documented MPMC group table, schema 1.

Uses typed integer GetValue access, never floating-point GetTuple for IDs.
This is application metadata in standard FieldData, not built-in VTK groups.
"""
from vtkmodules.vtkCommonCore import vtkUnsignedCharArray, vtkUnsignedIntArray, vtkUnsignedLongLongArray

NAMES = ('mpmc_group_schema_version', 'mpmc_group_location', 'mpmc_group_dimension',
         'mpmc_group_tag', 'mpmc_group_name_offsets', 'mpmc_group_name_utf8',
         'mpmc_group_member_offsets', 'mpmc_group_member_ids')
TYPES = (vtkUnsignedIntArray, vtkUnsignedCharArray, vtkUnsignedCharArray, vtkUnsignedIntArray,
         vtkUnsignedLongLongArray, vtkUnsignedCharArray, vtkUnsignedLongLongArray, vtkUnsignedLongLongArray)


def read_groups(grid):
    arrays = [grid.GetFieldData().GetArray(name) for name in NAMES]
    if all(array is None for array in arrays):
        return []
    values = []
    for name, array, kind in zip(NAMES, arrays, TYPES, strict=True):
        if array is None or not array.IsA(kind().GetClassName()) or array.GetNumberOfComponents() != 1:
            raise AssertionError('invalid VTK group array ' + name)
        values.append([int(array.GetValue(i)) for i in range(array.GetNumberOfTuples())])
    version, locations, dimensions, tags, name_offsets, name_bytes, member_offsets, members = values
    if version != [1] or any(len(data) != len(tags) for data in (locations, dimensions, name_offsets, member_offsets)):
        raise AssertionError('invalid group schema or row count')
    result = []
    nb = mb = 0
    for location, dimension, tag, ne, me in zip(locations, dimensions, tags, name_offsets, member_offsets, strict=True):
        if not (nb <= ne <= len(name_bytes) and mb <= me <= len(members)):
            raise AssertionError('invalid group offsets')
        result.append(dict(location=location, dimension=dimension, tag=tag,
                           name=bytes(name_bytes[nb:ne]).decode('utf-8'), members=members[mb:me]))
        nb, mb = ne, me
    if nb != len(name_bytes) or mb != len(members):
        raise AssertionError('unused group payload')
    return result


def physical_groups(grid):
    return {(group['dimension'], group['tag']): (group['name'], set(group['members']))
            for group in read_groups(grid)}


def write_groups(grid, groups):
    values = [[1], [], [], [], [], [], [], []]
    for group in groups:
        for index, key in ((1, 'location'), (2, 'dimension'), (3, 'tag')):
            values[index].append(group[key])
        values[5].extend(group['name'].encode('utf-8'))
        values[4].append(len(values[5]))
        values[7].extend(group['members'])
        values[6].append(len(values[7]))
    for name, kind, data in zip(NAMES, TYPES, values, strict=True):
        array = kind()
        array.SetName(name)
        for value in data:
            array.InsertNextValue(value)
        grid.GetFieldData().AddArray(array)
