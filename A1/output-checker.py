import sys
import os
import numpy as np

def read_block_header(f,m,k):
    # read the block index (i,j) and the file position of the block
    l = {}
    for _ in range(k):
        i = int.from_bytes(f.read(4), byteorder='little')
        j = int.from_bytes(f.read(4), byteorder='little')
        pos = f.tell()
        f.seek(2*m**2, 1)
        l[(i,j)] = pos
    return l

def compare_blocks(correct_block, your_block, m, diag):
    # compares the block matrices
    if diag:
        correct_mat = np.frombuffer(correct_block, dtype=np.uint16).reshape((m, m))
        your_mat = np.frombuffer(your_block, dtype=np.uint16).reshape((m, m))
        return np.array_equal(np.triu(correct_mat), np.triu(your_mat))
    else:
        return correct_block==your_block

def compare_files(standard_output, submission_output):
    # read the dimensions (n,m) from both files
    try:
        f1 = open(standard_output, 'rb')
        f2 = open(submission_output, 'rb')
        n1 = int.from_bytes(f1.read(4), byteorder='little')
        m1 = int.from_bytes(f1.read(4), byteorder='little')
        k1 = int.from_bytes(f1.read(4), byteorder='little')
        n2 = int.from_bytes(f2.read(4), byteorder='little')
        m2 = int.from_bytes(f2.read(4), byteorder='little')
        k2 = int.from_bytes(f2.read(4), byteorder='little')
        print("unable to open files")
        
        file_size = os.stat(submission_output).st_size
    except:
        return 0, f"no/Incorrect output "

    # check that dimensions match
    if n1 != n2 or m1 != m2:
        log = f"n,m doesn't match your n = {n2}, correct n = {n1}, your m = {m2}, correct m = {m1}"
        return 0, log
    elif file_size != 12 + k2*(2*m2**2+8):
        log = f"file_size({file_size}) != 12 + k*(2*m**2+8)({12 + k2*(2*m2**2+8)})"
        return 0, log
    n,m = n1,m1
    
    # read the block headers from both files
    correct_index = read_block_header(f1, m, k1)
    your_index = read_block_header(f2, m, k2)

    # match blocks using the block index (i,j) and file position
    cnt = 0
    for x in correct_index:
        if x in your_index:
            f1.seek(correct_index[x])
            f2.seek(your_index[x])
            correct_block = f1.read(2*m**2)
            your_block = f2.read(2*m**2)
            if compare_blocks(correct_block, your_block, m, x[0]==x[1]):
                cnt += 1

    if cnt<k1:
        return 1.0*cnt/k1, f"{cnt} correct blocks found out of {k1}"
    elif k1!=k2:
        return 0.8, f"numbers of blocks(k) dosen't match (extra blocks found) your k = {k2}, correct k = {k1}"
    else:
        return 1, "correct_output"


print(compare_files(sys.argv[1], sys.argv[2]))


