#!/usr/bin/python3

import os
import json
import sys
import argparse

script_root = "/home/student_id/training/Koryun/jacobi/PYTHON_VERSION"
class Error(Exception):
    pass

class NonValidRowError(Error):
    pass

class ExtraParametersError(Error):
    pass

class TooFewParametersError(Error):
    pass

class NoMatrixContentError(Error):
    pass

class Matrix:
	def __init__(self, matrix):
		self.matrix = matrix

	def print_data(self):
		print(self.matrix)

	def get_matrix(self):
		return self.matrix


"""
    the below function is responsible for getting filename from command line
"""
def arg_parse_foo():
    linear_parse=argparse.ArgumentParser(description="this script takes a file\
            as an argument from command line in which in ideal shold be a matrix\
            solves the matrix with Gaus method and \
            and afterwards generates another file containing the solution\
            of the given matrix")
    linear_parse.add_argument('-f', "--file", required = True)
    arguments = linear_parse.parse_args()
    return arguments.file


"""
    the below function is responsible for reading content from file input from command line
    and returning it as string
"""
def read_from_file(filename):
    try:
        with open(filename) as my_file:
            matrix_str = my_file.read()
            if 0 == len(matrix_str):
                try:
                    raise NoMatrixContentError
                except NoMatrixContentError:
                    print("file is empty !!\n")
                    sys.exit()
        return matrix_str
    except NoMatrixContentError:
        print("File not exist")
        sys.exit()


def generate_file(content, filename):
    with open(filename, 'w') as result_f:
        result_f.write(content)

def make_float_list(row_str_elem):
    for i in range(len(row_str_elem)):
        try:
            row_str_elem[i] = float(row_str_elem[i])
        except ValueError:
            print("not numeric content in file !!! \n")
            sys.exit()
    return row_str_elem


def define_raw_matrix(matrix_str):
    matrix = matrix_str.split("\n")
    matrix.pop()
    length = len(matrix)
    i = 0
    while i < length:
        matrix[i] = make_float_list(matrix[i].split())
        if len(matrix[i]) == matrix[i].count(0):
            matrix.remove(matrix[i])
            length -= 1 
            continue
        elif len(matrix[i]) - 1 == matrix[i].count(0) and 0 != matrix[i][-1]:
            try:
                raise NonValidRowError
            except NonValidRowError:
                print("the row is an invalid identity!\n")
                sys.exit()
        elif len(matrix[i]) != len(matrix[0]):
            try:
                raise NoMatrixContentError
            except:
                print("Invalid matrix or other content is in input file\n")
                sys.exit()
        i += 1
    return matrix

def is_Jacobi_method_appliable(A):
    if len(A[0]) != len(A) + 1:
        return True
    for i in range(len(A)):
        abs_sum = sum(map(abs, A[i])) - abs(A[i][i]) - abs(A[i][-1])
        if abs_sum > abs(A[i][i]):
            return True
    return False

def get_Jacobi_matrix(A):
    B = A
    for i in range(len(B)):
        x_i = B[i][i]
        del B[i][i]
        for j in range(len(B[i]) - 1):
            B[i][j] /= -x_i
        B[i][-1] /= x_i
    return B

def get_X_N_t_next(B, X_N_t):
    new_X_N_t = []
    tmp_list = list(X_N_t)
    for i in range(len(B)):
        tmp_list = list(X_N_t)
        new_X_N_t.append(float(0))
        del tmp_list[i]
        for j in range(len(tmp_list)):
            new_X_N_t[i] += tmp_list[j]*B[i][j]
        new_X_N_t[i] += B[i][-1]
    return new_X_N_t

def should_iteration_go_on(X_N_t, X_N_t_NEXT, aprox = 0.0001):
    for i in range(len(X_N_t)):
        if abs(abs(X_N_t_NEXT[i]) - abs(X_N_t[i])) > aprox:
            return True
    else:
        return False

def iterate_Jacobi(B, X_N_t):
    X_N_t_NEXT = get_X_N_t_next(B, X_N_t)
    while should_iteration_go_on(X_N_t, X_N_t_NEXT):
        tmp_list = list(X_N_t_NEXT)
        X_N_t_NEXT = get_X_N_t_next(B, X_N_t_NEXT)
        X_N_t = tmp_list
    return X_N_t_NEXT

def JACOBI_solve(matrix_obj):
    A = matrix_obj.get_matrix()
    if is_Jacobi_method_appliable(A):
        try:
            raise ValueError
        except ValueError:
            print("not possible to apply Jacobi method to given matrix\n")
            sys.exit()
    X_N_t = []
    if 1 == len(A):
        X_N_t = {}
        A[0][0] = A[0][1]/A[0][0]
        X_N_t["x_0"] = A[0][0]
        return X_N_t
    for i in range(len(A)):
        X_N_t.append(float(0))
    B = get_Jacobi_matrix(A)
    result_X_n = iterate_Jacobi(B, X_N_t)
    to_return = ""
    for elem in result_X_n:
        to_return += str(elem) + "   "
    return to_return

def remove_generated_files(filename):
    if os.path.exists(filename):
        os.remove(filename)
    else:
        print(f"The {filename} file does not exist")


def main():
    remove_generated_files(os.path.join(script_root, "/output.txt")) 
    remove_generated_files(os.path.join(script_root, "/result.txt")) 
    filename = arg_parse_foo()
    matrix_str = read_from_file(filename)
    raw_matrix = define_raw_matrix(matrix_str)
    matrix_obj = Matrix(raw_matrix)
    X_N_S = JACOBI_solve(matrix_obj)
    generate_file(X_N_S, "output.txt")

if __name__ == '__main__':
    main()
