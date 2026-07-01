import numpy as np


def to_shared_mut_matrix(parser):
    files = parser.files
    n_files = len(files)
    matrix = np.zeros((n_files, n_files), dtype=np.float64)

    for read_name in parser.reads_to_keep:
        read_data = parser.read_hash.get(read_name)
        if read_data is None:
            continue
        present_indices = [i for i, f in enumerate(files) if f in read_data]

        for ii in range(len(present_indices)):
            for jj in range(ii, len(present_indices)):
                i = present_indices[ii]
                j = present_indices[jj]
                matrix[i, j] += 1.0
                if i != j:
                    matrix[j, i] += 1.0

    return matrix


def to_matrix_presence_absence(parser):
    files = parser.files
    n_reads = len(parser.reads_to_keep)
    n_files = len(files)
    matrix = np.zeros((n_reads, n_files), dtype=np.int32)

    for ri, read_name in enumerate(parser.reads_to_keep):
        read_data = parser.read_hash.get(read_name)
        if read_data is None:
            continue
        for fi, f in enumerate(files):
            if f in read_data:
                matrix[ri, fi] = 1

    return matrix


def standardize_by_files_total(parser, matrix):
    total_reads_vector = parser.sample_total_reads_vector()
    if isinstance(matrix, np.ndarray):
        matrix = matrix.tolist()
    n_cols = len(matrix[0]) if matrix else 0

    result = np.zeros((len(matrix), n_cols), dtype=np.float64)
    for i, row in enumerate(matrix):
        if total_reads_vector[i] == 0:
            result[i, :] = 0.0
        else:
            for j in range(n_cols):
                if total_reads_vector[j] > 0:
                    result[i, j] = row[j] / total_reads_vector[j]
                else:
                    result[i, j] = 0.0
    return result


def standardize_by_files_total_assymetric(parser, matrix):
    total_reads_vector = parser.sample_total_reads_vector()
    n_cols = len(matrix[0]) if matrix else 0

    result = np.zeros((len(matrix), n_cols), dtype=np.float64)
    for i, row in enumerate(matrix):
        if total_reads_vector[i] == 0:
            result[i, :] = 0.0
        else:
            for j in range(n_cols):
                result[i, j] = row[j] / total_reads_vector[i]
    return result


def invert_matrix(matrix):
    if isinstance(matrix, np.ndarray):
        return 1.0 - matrix
    return [[1.0 - val for val in row] for row in matrix]


def transpose_matrix(matrix):
    if isinstance(matrix, np.ndarray):
        return matrix.T
    if not matrix:
        return []
    return list(map(list, zip(*matrix)))
