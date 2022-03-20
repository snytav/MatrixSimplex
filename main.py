# https://towardsdatascience.com/developing-the-simplex-method-with-numpy-and-matrix-operations-16321fd82c85
# https://in.mathworks.com/matlabcentral/fileexchange/85223-linear-programming-simplex-algorithm
# https://gist.github.com/snytav/b6d59126610325a0b74562d09a343f3d
import numpy as np

def init_simplex(A):
    # sizes of basic and nonbasic vectors
    basicSize = A.shape[0]  # number of constraints, m
    nonbasicSize = A.shape[1] - basicSize  # n-m, number of variables

    # global index tracker of variables of basic and nonbasic variables (objective)
    # that is, index 0 corresponds with x_0, 1 with x_1 and so on.  So each index corresponds with a variable
    cindx = [i for i in range(0, len(c))]

    # basic variable coefficients
    cbT = np.array(c[nonbasicSize:])

    # nonbasic variable coefficients
    cnT = np.array(c[:nonbasicSize])
    return cbT,cnT,cindx,nonbasicSize

def get_non_basic_ratios(Ahat,bHat):
    # now we want to iterate through Ahat and bHat and pick the minimum ratios
    # only take ratios of variables with Ahat_i values greater than 0
    # pick smallest ratio to get variable that will become nonbasic.
    ratios = []
    for i in range(0, len(bHat)):
        Aval = Ahat[i]
        Bval = bHat[i]

        # don't look at ratios with val less then or eqaul to 0, append to keep index
        if (Aval <= 0):
            ratios.append(10000000)
            continue
        ratios.append(Bval / Aval)

    return ratios

def switch_indices(ratios,cnT,cbT,cindx,nonbasicSize,cnMinIndx):
    ratioMinIndx = np.argmin(ratios)

    # switch basic and nonbasic variables using the indices.
    cnT[cnMinIndx], cbT[ratioMinIndx] = cbT[ratioMinIndx], cnT[cnMinIndx]
    # switch global index tracker indices
    cindx[cnMinIndx], cindx[ratioMinIndx + nonbasicSize] = cindx[ratioMinIndx + nonbasicSize], cindx[cnMinIndx]
    return cnT,cbT,cindx

def core_simplex(cindx,nonbasicSize,A,cnT,cbT):
    # run core simplex method until reach the optimal solution
    while True:

        # keep track of current indices of basic and non-basic variables
        cbIndx = cindx[nonbasicSize:]
        cnIndx = cindx[:nonbasicSize]

        # basis matrix
        B = A[:, cbIndx]
        Binv = np.linalg.inv(B)

        # nonbasic variable matrix
        N = A[:, cnIndx]

        # bHat, the values of the basic variables
        # recall that at the start the basic variables are the slack variables, and
        # have values equal the vector b (as primary variables are set to 0 at the start)
        bHat = Binv @ b
        yT = cbT @ Binv

        # use to check for optimality, determine variable to enter basis
        cnHat = cnT - (yT @ N)

        # find indx of minimum value of cnhat, this is the variable to enter the basis
        cnMinIndx = np.argmin(cnHat)

        # break out of loop, returning values if all values of cnhat are above 0
        if (all(i >= 0 for i in cnHat)):
            # use cbIndx to get index values of variables in bHat, and the corresponding index
            # values in bHat are the final solution values for each of the corresponding variables
            # ie value 0 in dbIndx corresponds with first variable, so whatever the index for the 0 is
            # is the index in bHat that has the solution value for that variable.
            return cbT, cbIndx, cnT, cnIndx, bHat, cnHat

        # this is the index for the column of coeffs in a for the given variable
        indx = cindx[cnMinIndx]

        Ahat = Binv @ A[:, indx]

        ratios = get_non_basic_ratios(Ahat, bHat)
        # now we want to iterate through Ahat and bHat and pick the minimum ratios
        # only take ratios of variables with Ahat_i values greater than 0
        # pick smallest ratio to get variable that will become nonbasic.
        # ratios = []
        # for i in range(0, len(bHat)):
        #     Aval = Ahat[i]
        #     Bval = bHat[i]
        #
        #     # don't look at ratios with val less then or eqaul to 0, append to keep index
        #     if (Aval <= 0):
        #         ratios.append(10000000)
        #         continue
        #     ratios.append(Bval / Aval)

        cnT,cbT,cindx = switch_indices(ratios, cnT, cbT, cindx, nonbasicSize, cnMinIndx)
        # ratioMinIndx = np.argmin(ratios)
        #
        # # switch basic and nonbasic variables using the indices.
        # cnT[cnMinIndx], cbT[ratioMinIndx] = cbT[ratioMinIndx], cnT[cnMinIndx]
        # # switch global index tracker indices
        # cindx[cnMinIndx], cindx[ratioMinIndx + nonbasicSize] = cindx[ratioMinIndx + nonbasicSize], cindx[cnMinIndx]
        # # now repeat the loop



def Simplex(A, b, c):
    '''Takes input vars, computs corresponding values,
    then uses while loop to iterate until a basic optimal solution is reached.
    RETURNS: cbT, cbIndx, cnT, cnIndx, bHat, cnHat.
    cbT, cbIndex is final basic variable values, and indices
    cnT, cnIndex is final nonbasic variable values and indices
    bHat is final solution values,
    cnHat is optimality condition'''

    # # sizes of basic and nonbasic vectors
    # basicSize = A.shape[0]  # number of constraints, m
    # nonbasicSize = A.shape[1] - basicSize  # n-m, number of variables
    #
    # # global index tracker of variables of basic and nonbasic variables (objective)
    # # that is, index 0 corresponds with x_0, 1 with x_1 and so on.  So each index corresponds with a variable
    # cindx = [i for i in range(0, len(c))]
    #
    # # basic variable coefficients
    # cbT = np.array(c[nonbasicSize:])
    #
    # # nonbasic variable coefficients
    # cnT = np.array(c[:nonbasicSize])
    cbT, cnT, cindx, nonbasicSize = init_simplex(A)
    # run core simplex method until reach the optimal solution
    return core_simplex( cindx, nonbasicSize, A, cnT, cbT)
    # while True:
    #
    #     # keep track of current indices of basic and non-basic variables
    #     cbIndx = cindx[nonbasicSize:]
    #     cnIndx = cindx[:nonbasicSize]
    #
    #     # basis matrix
    #     B = A[:, cbIndx]
    #     Binv = np.linalg.inv(B)
    #
    #     # nonbasic variable matrix
    #     N = A[:, cnIndx]
    #
    #     # bHat, the values of the basic variables
    #     # recall that at the start the basic variables are the slack variables, and
    #     # have values equal the vector b (as primary variables are set to 0 at the start)
    #     bHat = Binv @ b
    #     yT = cbT @ Binv
    #
    #     # use to check for optimality, determine variable to enter basis
    #     cnHat = cnT - (yT @ N)
    #
    #     # find indx of minimum value of cnhat, this is the variable to enter the basis
    #     cnMinIndx = np.argmin(cnHat)
    #
    #     # break out of loop, returning values if all values of cnhat are above 0
    #     if (all(i >= 0 for i in cnHat)):
    #         # use cbIndx to get index values of variables in bHat, and the corresponding index
    #         # values in bHat are the final solution values for each of the corresponding variables
    #         # ie value 0 in dbIndx corresponds with first variable, so whatever the index for the 0 is
    #         # is the index in bHat that has the solution value for that variable.
    #         return cbT, cbIndx, cnT, cnIndx, bHat, cnHat
    #
    #     # this is the index for the column of coeffs in a for the given variable
    #     indx = cindx[cnMinIndx]
    #
    #     Ahat = Binv @ A[:, indx]
    #
    #     ratios = get_non_basic_ratios(Ahat,bHat)
    #     # now we want to iterate through Ahat and bHat and pick the minimum ratios
    #     # only take ratios of variables with Ahat_i values greater than 0
    #     # pick smallest ratio to get variable that will become nonbasic.
    #     # ratios = []
    #     # for i in range(0, len(bHat)):
    #     #     Aval = Ahat[i]
    #     #     Bval = bHat[i]
    #     #
    #     #     # don't look at ratios with val less then or eqaul to 0, append to keep index
    #     #     if (Aval <= 0):
    #     #         ratios.append(10000000)
    #     #         continue
    #     #     ratios.append(Bval / Aval)
    #
    #     ratioMinIndx = np.argmin(ratios)
    #
    #     # switch basic and nonbasic variables using the indices.
    #     cnT[cnMinIndx], cbT[ratioMinIndx] = cbT[ratioMinIndx], cnT[cnMinIndx]
    #     # switch global index tracker indices
    #     cindx[cnMinIndx], cindx[ratioMinIndx + nonbasicSize] = cindx[ratioMinIndx + nonbasicSize], cindx[cnMinIndx]
    #     # now repeat the loop

# Example Input
A = np.array([[-2, 1, 1, 0, 0],
             [-1, 2, 0, 1, 0],
             [1, 0, 0, 0, 1]])

b = np.array([2, 7, 3])

c = np.array([-1, -2, 0, 0, 0])



cbT, cbIndx, cnT, cnIndx, bHat, cnHat = Simplex(A, b, c)
print(cbT)
print(cbIndx)
print(cnT)
print(cnIndx)
print(bHat)
print(cnHat)
qq = 0