######################################################################
# Name: Jose L. Camacho
# Class: CS 3010
# Professor = Edwin Rodríguez
# Description - Proof of concept for naive gaussean elemination
# for a system of equations.
######################################################################
import sys  # needed for command line arguements


# Naive Gaussian Elimination
def FwdElimination(coeff: list, const: list):
    # since coeff is essentially an nxn matrix to make it easier on me and the reader
    # i've assigned n the size of the matrix
    n = len(coeff)

    for k in range(n - 1):
        for i in range(k + 1, n):
            mult = coeff[i][k] / coeff[k][k]
            for j in range(k + 1, n):
                coeff[i][j] = coeff[i][j] - mult * coeff[k][j]
            const[i] = const[i] - mult * const[k]
    return coeff, const


def BackSubstitution(coeff: list, const: list, sol: list):
    n = len(coeff)
    sol[-1] = const[-1] / coeff[-1][-1]
    for i in range(n - 1, -1, -1):
        sum = const[i]
        for j in range(i + 1, n):
            sum = sum - coeff[i][j] * sol[j]
        sol[i] = sum / coeff[i][i]
    return sol


def NaiveGaussian(coeff: list, const: list) -> list:
    sol = [0 for x in range(len(coeff))]

    coeff, cosnt = FwdElimination(coeff, const)
    sol = BackSubstitution(coeff, const, sol)
    return sol


# SPP Gaussian Algorithm
def SPPFwdElimination(coeff: list, const: list, ind: list):
    # list containing scaling factors
    scaling = [0 for x in range(len(coeff))]
    n = len(coeff)

    for i in range(n):
        smax = 0
        for j in range(n):
            smax = max(smax, abs(coeff[i][j]))
        scaling[i] = smax

    for k in range(n - 1):
        rmax = 0
        maxind = k

        for i in range(k, n):
            # ratio of coefficient to scaling factor
            r = abs(coeff[ind[i]][k] / scaling[ind[i]])
            if r > rmax:
                rmax = r
                maxind = i
        ind[maxind], ind[k] = ind[k], ind[maxind]  # swap values

        for i in range(k + 1, n):
            mult = coeff[ind[i]][k] / coeff[ind[k]][k]
            for j in range(k + 1, n):
                coeff[ind[i]][j] = coeff[ind[i]][j] - mult * coeff[ind[k]][j]
            const[ind[i]] = const[ind[i]] - mult * const[ind[k]]

    return coeff, const, ind


def SPPBackSubst(coeff: list, const: list, sol: list, ind: list):
    n = len(coeff)
    sol[-1] = const[ind[-1]] / coeff[ind[-1]][-1]

    for i in range(n - 1, -1, -1):
        sum = const[ind[i]]
        for j in range(i + 1, n):
            sum = sum - coeff[ind[i]][j] * sol[j]
        sol[i] = sum / coeff[ind[i]][i]
    return sol


def SPPGaussian(coeff: list, const: list) -> list:
    sol = [0 for x in range(len(coeff))]
    ind = [0 for x in range(len(coeff))]
    n = len(coeff)
    for i in range(n):
        ind[i] = i

    coeff, const, ind = SPPFwdElimination(coeff, const, ind)
    sol = SPPBackSubst(coeff, const, sol, ind)
    return sol


def cline(buff: list) -> list:
    for i in range(len(buff)):
        buff[i] = float(buff[i])
    return buff


def clean(buff: list) -> list:
    nbuff = list()
    for i in range(len(buff)):
        if buff[i] != '':
            nbuff.append(buff[i])
    return nbuff


def main() -> None:
    flag = 0  # flag to check if --spp was added
    coeff = []  # coefficient matrix  (a)
    const = list()  # constant matrix (b) in Ax = b
    filename = ''  # filename
    q = 0  # counter for file processing

    argList = sys.argv  # assign cmd line args to dedicated var

    # parse arguments list
    for i in range(1, len(argList)):
        if argList[i] == '--spp':
            flag = 1
        if '.lin' in argList[i]:
            filename = argList[i]

    # begin processing input file
    if filename == '':
        print('ERROR no input file')
        sys.exit()
    else:
        with open(filename, 'r', encoding='utf-8') as f:
            # assume first line is the dimension of the matrix and final line is the values
            # the matrix is equal to; i.e. Ax = b
            for line in f:
                line = line.strip('\n')
                s = line.split(' ')
                s = clean(s)
                if q == 0:
                    n = int(s[0])
                    q += 1
                elif q <= n + 1:
                    coeff.append(cline(s))
                else:
                    const = cline(s)
                q += 1

    if flag:
        sol = SPPGaussian(coeff, const)
        filename = 'spp_' + filename
    else:
        sol = NaiveGaussian(coeff, const)

    # modify filename to have '.sol' extension
    filename = filename.replace('.lin', '.sol')

    # convert sol to string
    sol = ' '.join([str(i) for i in sol])
    sol += '\n'

    # generate .sol file
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(sol)
    except IOError:
        print('I/O error encountered')


# program execution begins here
if __name__ == '__main__':
    main()
