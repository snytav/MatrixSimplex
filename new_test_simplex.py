import numpy as np
from torch_simplex import A,b,c,Simplex










if __name__ == '__main__':
    import numpy as np
    # a = np.array([0, 1, 2])
    # print(np.tile(a, 2))
    # print(np.tile(a, (2, 2)))
    AA = np.tile(A,(2,2))
    bb = np.tile(b,(1,2))
    cc = np.tile(c,(2,1))


    cbT, cbIndx, cnT, cnIndx, bHat, cnHat = Simplex(A, b, c)
    print(cbT)
    print(cbIndx)
    print(cnT)
    print(cnIndx)
    print(bHat)
    print(cnHat)