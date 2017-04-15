################################################################################
#
# Author: Christian Jost (jost@hiskp.uni-bonn.de)
# Date:   Februar 2015
#
# Copyright (C) 2015 Christian Jost
# 
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with tmLQCD. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
#
# Function: Implements functions for manipulating correlation function matrices.
#
# For informations on input parameters see the description of the function.
#
################################################################################

__all__ = ["create_corr_matrix_cosh", "create_corr_matrix_sinh", "calculate_gevp", "shift_corr_matrix", "write_eigenvals_cosh", "write_eigenvals_sinh", "reorder_by_evec", "reorder_by_eval", "reorder_by_evec_all", "reorder_by_evec_all_norm"]

import os
import numpy as np
import scipy.linalg as spla
import input_output as io
import bootstrap
import itertools

def shift_corr_matrix(cmat, dt, dE=None, axis=1):
    """Weight-shift the correlation function matrix.

    This is based on the paper by Dudek et al, Phys.Rev. D86, 034031 (2012).
    First the matrix is weighted by exp(dE*t) on every timeslice and then shifted.
    If dE is not given, the matrix is only shifted

    Args:
        cmat: correlation function matrix
        dt: number of timeslices to shift
        dE: the energy for the weighting
        axis: axis along which the shift is performed

    Returns:
        shifted correlation function matrix
    """
    # check whether axis is in the number of dimension
    ndim = cmat.ndim
    if axis > ndim or axis < 0:
        print("axis not in range of the dimensions. Not doing anything")
        return cmat

    # calculate shape of the new matrix
    cmshape = np.asarray(cmat.shape)

    # weighting of the matrix
    if dE:
        for t in xrange(cmshape[axis]):
          # TODO(CJ): Think of a way to do this more general
          cmat[:,t] = np.exp(dE*t) * cmat[:,t]

    # if dt is zero, don't shift
    if dt is 0:
        return cmat

    # create the new array
    cm = np.zeros(cmshape)

    # fill the new array
    for i in xrange(cmshape[axis]):
        if i+dt >= cmshape[axis]:
          cm[:,i] = cmat[:,i] - cmat[:,(i+dt)- cmshape[axis]]
        else:
          cm[:,i] = cmat[:,i] - cmat[:,i+dt]


    # return shifted matrix
    return cm

def create_corr_matrix_cosh(matrix_index, filepath, fileprefix, filesuffix,
                       complex=True, verbose=0):
    """Creates a correlation function matrix.

    Reads different correlation functions and inserts them into a matrix. The
    matrix is filled row majored, the correlation functions matrix is stored
    column majored. It is assumed that the matrix is symmetric, the off
    diagonal elements are symmetrized.

    Args:
        nbsamples: Number of bootstrap samples created.
        filepath: The path to the data, including the fileprefix.
        filestring: A list of the changing parts of the filenames. The length
                    of the list gives the size of the matrix.
        filesuffix: The suffix of the data files.
        column: The column of the input file to be read. The same column is
                read from every file!
        verbose: Changes the amount of information printed.

    Returns:
        A numpy array with four axis. The first axis is the column of operators,
        the second axis is the row of operators, the third axis is the number
        of the bootstrap sample and the last axis is the time dependence.
        The second returned argument is the time extend of the correlation
        function matrix, that is (T/2)+1, where T is the length of the original
        correlation function in time.
    """
 
    fname = "%s/%s%i%i%s" % (filepath, fileprefix, 0, 0, filesuffix)
    datashape = io.read_header(fname)
    nsample = datashape[0]
    Lt = datashape[1]
    type = datashape[2]
    if type == 1:
       data_matrix = np.zeros((nsample, Lt, len(matrix_index), len(matrix_index)), dtype="complex")
    else:
       data_matrix = np.zeros((nsample, Lt, len(matrix_index), len(matrix_index)), dtype="float")
      
    for i in range(len(matrix_index)):
        for j in range(len(matrix_index)):
          fname = "%s/%s%i%i%s" % (filepath, fileprefix, matrix_index[i], matrix_index[j], filesuffix)
          data = io.read_data_ascii(fname)
          if data.shape[0] != nsample or data.shape[1] != Lt:
             print("incompatible correlators")
          else:
             data_sys = bootstrap.sym(data)
             data_matrix[:,:,i,j] = data_sys

    corr_mat_symm = np.zeros_like(data_matrix)
    for _s in range(0, nsample):
        for _t in range(0, Lt):
            corr_mat_symm[_s, _t] = (data_matrix[_s, _t] + data_matrix[_s, _t].T) / 2.

    if not complex and type == 1:
       return np.real(corr_mat_symm)
    else:
       return corr_mat_symm

def create_corr_matrix_sinh(matrix_index, filepath, fileprefix, filesuffix,
                       complex=True, verbose=0):
    """Creates a correlation function matrix.

    Reads different correlation functions and inserts them into a matrix. The
    matrix is filled row majored, the correlation functions matrix is stored
    column majored. It is assumed that the matrix is symmetric, the off
    diagonal elements are symmetrized.

    Args:
        nbsamples: Number of bootstrap samples created.
        filepath: The path to the data, including the fileprefix.
        filestring: A list of the changing parts of the filenames. The length
                    of the list gives the size of the matrix.
        filesuffix: The suffix of the data files.
        column: The column of the input file to be read. The same column is
                read from every file!
        verbose: Changes the amount of information printed.

    Returns:
        A numpy array with four axis. The first axis is the column of operators,
        the second axis is the row of operators, the third axis is the number
        of the bootstrap sample and the last axis is the time dependence.
        The second returned argument is the time extend of the correlation
        function matrix, that is (T/2)+1, where T is the length of the original
        correlation function in time.
    """
 
    fname = "%s/%s%i%i%s" % (filepath, fileprefix, 0, 0, filesuffix)
    datashape = io.read_header(fname)
    nsample = datashape[0]
    Lt = datashape[1]
    type = datashape[2]
    if type == 1:
       data_matrix = np.zeros((nsample, Lt, len(matrix_index), len(matrix_index)), dtype="complex")
    else:
       data_matrix = np.zeros((nsample, Lt, len(matrix_index), len(matrix_index)), dtype="float")
      
    for i in range(len(matrix_index)):
        for j in range(len(matrix_index)):
          fname = "%s/%s%i%i%s" % (filepath, fileprefix, matrix_index[i], matrix_index[j], filesuffix)
          data = io.read_data_ascii(fname)
          if data.shape[0] != nsample or data.shape[1] != Lt:
             print("incompatible correlators")
          else:
             data_matrix[:,:,i,j] = data

    corr_mat_symm = np.zeros_like(data_matrix)
    for _s in range(0, nsample):
        for _t in range(0, Lt):
            corr_mat_symm[_s, _t] = (data_matrix[_s, _t] + data_matrix[_s, _t].T) / 2.

    if not complex and type == 1:
       return np.real(corr_mat_symm)
    else:
       return corr_mat_symm

      
def reorder_by_eval(eigen_values, eigen_vectors, t0):
    nsamples = eigen_values.shape[0]
    T = eigen_values.shape[1]
    dim = eigen_values.shape[2]

    order = np.zeros((T,dim), dtype="int")
    for _n in range(nsamples):
       for _t in range(T):
          if _t >= t0:
            order[_t] = np.array(list(reversed(sorted(range(dim), key = eigen_values[_n, _t].__getitem__))))
          else:
            order[_t] = np.array(list(sorted(range(dim), key = eigen_values[_n, _t].__getitem__)))
          eigen_values[_n, _t] = eigen_values[_n, _t][order[_t]]
          eigen_vectors[_n, _t] = eigen_vectors[_n, _t][order[_t]]
    return 


def reorder_by_evec(eigen_values, eigen_vectors, Ct0, t0):
    print(eigen_values.shape)
    print(eigen_vectors.shape)
    nsamples = eigen_values.shape[0]
    T = eigen_values.shape[1]
    dim = eigen_values.shape[2]
    print(Ct0.shape)

    order = np.zeros((T,dim), dtype="int")
#    for _n in range(nsamples):
    for _n in range(1):
#    for _n in range(1):
# sort t0+1 and t0-1 by eigen values
       order[t0+1] = np.array(list(reversed(sorted(range(dim), key = eigen_values[_n, t0+1].__getitem__))))

       for _t in range(t0+2, T):
         test1 = np.dot(Ct0[_n], eigen_vectors[_n, _t])
         test2 = np.dot(eigen_vectors[_n, _t].transpose(), test1)
         print(test2)
         tref = t0+1
         for _i2 in range(dim):
             Ct0ev2 = np.dot(Ct0[_n], eigen_vectors[_n, _t, :, _i2])
             ev1Ct0ev2 = np.abs(np.dot(eigen_vectors[_n, tref, :, range(dim)], Ct0ev2))
#             ev1Ct0ev2_check = np.abs(np.dot(eigen_vectors[_n, _t, :, 0], Ct0ev2))
             maxi = np.where(ev1Ct0ev2 == ev1Ct0ev2.max())[0][0]
             order[_t, _i2] = order[tref, maxi]
             print("------------sample: %i, t: %i, eig%i------------------" % (_n, _t, _i2))
             print(Ct0ev2)
             print(ev1Ct0ev2)
#             print(ev1Ct0ev2_check)
         print(order[_t])
         while sorted(order[_t]) != range(dim):
           tref = tref+1
           if tref >= _t:
             print("sample: %i, t: %i sort by eigen vector failed, try sorting by eigen values..." % (_n, _t)) 
             order[_t] = np.array(list(reversed(sorted(range(dim), key = eigen_values[_n, _t].__getitem__))))         
#           print("sample:%i, t: %i, choosing other tref..." % (_n, _t))
           else:
             for _i2 in range(dim):
               Ct0ev2 = np.dot(Ct0[_n], eigen_vectors[_n, _t, :, _i2])
               ev1Ct0ev2 = np.abs(np.dot(eigen_vectors[_n, tref, :, range(dim)], Ct0ev2[0]))
               maxi = np.where(ev1Ct0ev2 == ev1Ct0ev2.max())[0][0]
               order[_t, _i2] = order[tref, maxi]
 
       for _t in range(t0, 0, -1):
           tref = t0+1
           for _i2 in range(dim):
             Ct0ev2 = np.dot(Ct0[_n], eigen_vectors[_n, _t, :,  _i2])
             ev1Ct0ev2 = np.abs(np.dot(eigen_vectors[_n, tref, :, range(dim)], Ct0ev2[0]))
             maxi = np.where(ev1Ct0ev2 == ev1Ct0ev2.max())[0][0]
             order[_t, _i2] = order[_t+1, maxi]
           for _check in range(dim):
             if _check not in order[_t]:
               print("sample: %i, t: %i sort by eigen vector failed, try sorting by eigen values..." % (_n, _t)) 
               order[_t] = np.array(list(sorted(range(dim), key = eigen_values[_n, _t].__getitem__))) 

       for _t in range(1, T):  
         eigen_values[_n, _t] = eigen_values[_n, _t][order[_t]]
         eigen_vectors[_n, _t] = eigen_vectors[_n, _t][order[_t]]
             
    return

def reorder_by_evec_all(eigen_values, eigen_vectors, Ct0, t0):
    nsamples = eigen_values.shape[0]
    T = eigen_values.shape[1]
    dim = eigen_values.shape[2]

    order = np.zeros((T,dim), dtype="int")
    for _n in range(nsamples):
       order[t0+1] = np.array(list(reversed(sorted(range(dim), key = eigen_values[_n, t0+1].__getitem__))))

       for _t in range(T):
         if _t == t0+1:
            continue

         det_max = 0.
         for _perm in itertools.permutations(range(dim)):
            evec = eigen_vectors[_n, t0+1]
            evec_product = np.identity(dim)
            for _i in range(dim):
              evec[:, _perm[_i]] = eigen_vectors[_n, _t, :, _i]
              evec_product = np.dot(evec_product, evec)
            det = np.abs(np.linalg.det(evec_product))
            if det > det_max:
               det_max = det
               order[_t] = order[t0+1][_perm]

       for _t in range(1, T):  
         eigen_values[_n, _t] = eigen_values[_n, _t][order[_t]]
         eigen_vectors[_n, _t] = eigen_vectors[_n, _t][order[_t]]
             
    return

       
def reorder_by_evec_all_norm(eigen_values, eigen_vectors, Ct0, t0):
    nsamples = eigen_values.shape[0]
    T = eigen_values.shape[1]
    dim = eigen_values.shape[2]

    Ct02 = np.zeros_like(Ct0)
    eigen_vectors_norm = np.zeros_like(eigen_vectors)
    for _n in range(nsamples):
      Ct02[_n] = spla.sqrtm(Ct0[_n])
      for _t in range(T):
        eigen_vectors_norm[_n, _t] = np.dot(Ct02[_n], eigen_vectors[_n, _t])

    
    order = np.zeros((T,dim), dtype="int")
    for _n in range(nsamples):
       order[t0+1] = np.array(list(reversed(sorted(range(dim), key = eigen_values[_n, t0+1].__getitem__))))

       for _t in range(T):
         if _t == t0+1:
            continue

         det_max = 0.
         for _perm in itertools.permutations(range(dim)):
            evec = eigen_vectors_norm[_n, t0+1]
            evec_product = np.identity(dim)
            for _i in range(dim):
              evec[:, _perm[_i]] = eigen_vectors_norm[_n, _t, :, _i]
              evec_product = np.dot(evec_product, evec)
            det = np.abs(np.linalg.det(evec_product))
            if det > det_max:
               det_max = det
               order[_t] = order[t0+1][_perm]

       for _t in range(1, T):  
         eigen_values[_n, _t] = eigen_values[_n, _t][order[_t]]
         eigen_vectors[_n, _t] = eigen_vectors[_n, _t][order[_t]]
             
    return



 
def solve_gevp_gen(a, t_0, sym=True):
    """Generator that returns the eigenvalues for t_0 -> t where t is in
       (t_0, t_max].
       
       Calculate the eigenvalues of the generalised eigenvalue problem using
       the scipy.linalg.eigh solver.
    
       Args:
           a: The matrix that is used.
           t_0: The timeslice of the inverted matrix.

       Returns:
           A list of the eigenvalues, a list of the eigenvectors and the time
           slice on which the eigensystem was calculated.
    """
    # B is the matrix at t=t_0
    B = a[t_0]
    dim = B.shape[0]
    if sym:
       T2 = a.shape[0]/2 + 1;
    else:
       T2 = a.shape[0]/2;

    # define the eigensystem solver function as a lambda function
    try:
        f = lambda A: spla.eigh(a=A, b=B)
    except LinAlgError:
        return

    # initialization
    count = 0
    # calculate the eigensystem for t in (t_0, T/2+1)
    for j in range(0, T2):
#        try:
            # calculate the eigensystems
            eigenvalues, eigenvectors = f(a[j])
            count += 1
            yield eigenvalues, eigenvectors, j

#        except (spla.LinAlgError, TypeError) as e:
#            print("solve eigen failed")
#            return

def calculate_gevp(m, t0=1, sym=True):
    """Solves the generalized eigenvalue problem of a correlation function
    matrix.

    The function takes a bootstrapped correlation function matrix and calculates
    the eigenvectors and eigenvalues of the matrix. The algorithm relies on the
    matrix being symmetric or hermitian. The matrix m should have 4 axis, as
    laid out in corr_matrix.py.
    The eigenvectors are calculated but not stored.

    Args:
        m: The data in a numpy array.
        t0: The timeslice used for the inversion.

    Returns:
        A numpy array with three axis. The first axis is the bootstrap sample
        number, the second axis is the time, the third axis is the eigenvalue
        numbers. The time extend is the same as in original data, but the times
        up to t0 are filled with zeros.
    """
    # Initialize the eigenvalue array
    if sym:
       T2 = m.shape[1]/2 + 1
    else:
       T2 = m.shape[1]/2
    values_array = np.zeros((m.shape[0], T2, m.shape[2]))
    vectors_array = np.zeros((m.shape[0], T2, m.shape[2], m.shape[2]), dtype="complex")
    # iterate over the bootstrap samples
    for _samples in range(0, m.shape[0]):
#    for _samples in range(1):
        # iterate over the eigensystems
        for eigenvalues, eigenvectors, _t in solve_gevp_gen(m[_samples], t0, sym):
            # save the eigenvalues to the array
            values_array[_samples, _t] = eigenvalues
            vectors_array[_samples, _t] = eigenvectors
    return values_array, vectors_array

def write_eigenvals_cosh(fname, eigvals):
   """wirte the eigvals to a text file in the same format as corr
      Args:
         eigvals: three dimensional array. the first axis is the sample number, the second is time,          the third axis is the number of eigenvalues. 
   """

   nsample = eigvals.shape[0]
   T2 = eigvals.shape[1]
   nvals = eigvals.shape[2]
   if len(fname) != nvals:
      print("number of output filenames and number of eigen values do not match.")
      os.sys.exit(-1)
   Lt = (T2-1)*2
   _eigoutput = np.zeros((nsample, Lt, nvals))
   for _t in range(T2):
         _eigoutput[:,_t,:] = eigvals[:,_t,:]
   for _t in range(T2, Lt):
         _eigoutput[:,_t,:] = eigvals[:,Lt-_t,:]
   for i in range(nvals):
     io.write_data_ascii(_eigoutput[:,:,i], fname[i], complex=False)        

       
def write_eigenvals_sinh(fname, eigvals):
   """wirte the eigvals to a text file in the same format as corr
      Args:
         eigvals: three dimensional array. the first axis is the sample number, the second is time,          the third axis is the number of eigenvalues. 
   """

   nsample = eigvals.shape[0]
   T2 = eigvals.shape[1]
   nvals = eigvals.shape[2]
   if len(fname) != nvals:
      print("number of output filenames and number of eigen values do not match.")
      os.sys.exit(-1)
   Lt = T2*2
   _eigoutput = np.zeros((nsample, Lt, nvals))
   for _t in range(T2):
         _eigoutput[:,_t,:] = eigvals[:,_t,:]
   for _t in range(T2, Lt):
         _eigoutput[:,_t,:] = -1.0*eigvals[:,Lt-_t-1,:]
   for i in range(nvals):
     io.write_data_ascii(_eigoutput[:,:,i], fname[i], complex=False)        


