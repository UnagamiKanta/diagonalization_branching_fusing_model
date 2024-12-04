import numpy as np
from scipy.linalg import eig

def calc_eigen_value_for_two_site(alpha, beta, gamma, delta, D):

    matrix = np.array([
        [0,           -delta,           -delta,             -2*alpha],
        [0, beta + gamma + D,               -D,               -gamma],
        [0,               -D, beta + gamma + D,               -gamma],
        [0,            -beta,            -beta, 2 * alpha + 2 *gamma]
    ])

    eigen_values = eig(matrix)[0]
    
    return eigen_values

def calc_eigen_value_for_three_site(alpha, beta, gamma, delta, D):
    A = 2 * alpha


    matrix = np.array([
        [0,                   -delta,                   -delta,                   -delta,                                            -A,                                           -A,                                            -A,                          0],
        [0, 2 * beta + delta + 2 * D,                       -D,                       -D,                                -gamma - delta,                               -gamma - delta,                                             0,                         -A],
        [0,                       -D, 2 * beta + delta + 2 * D,                       -D,                                -gamma - delta,                                            0,                                -gamma - delta,                         -A],
        [0,                       -D,                       -D, 2 * beta + delta + 2 * D,                                             0,                               -gamma - delta,                                -gamma - delta,                         -A],
        [0,                    -beta,                    -beta,                        0,  A + 2 * beta + 2 * gamma + 2 * delta + 2 * D,                                           -D,                                            -D,                     -gamma],
        [0,                    -beta,                        0,                    -beta,                                            -D, A + 2 * beta + 2 * gamma + 2 * delta + 2 * D,                                            -D,                     -gamma],
        [0,                        0,                    -beta,                    -beta,                                            -D,                                           -D,  A + 2 * beta + 2 * gamma + 2 * delta + 2 * D,                     -gamma],
        [0,                        0,                        0,                        0,                                     -2 * beta,                                    -2 * beta,                                     -2 * beta,      6 * alpha + 3 * gamma]
    ])


    eigen_values = eig(matrix)[0]
    
    return eigen_values


def are_all_real_numbers(eigen_values):
    return all([np.isclose(c.imag, 0, atol=tolerance) for c in eigen_values])

if __name__ == '__main__':
    tolerance = 1e-15  # 許容誤差の設定
    for alpha in np.arange(0, 1, 0.1):
        for beta in np.arange(0, 1, 0.1):
            for gamma in np.arange(0, 1, 0.1):
                for delta in np.arange(0, 1, 0.1):
                    for D in np.arange(0, 1, 0.1):
                        if alpha + beta + gamma + delta + D > 0.5:
                            continue
                        eigen_values = calc_eigen_value_for_three_site(alpha, beta, gamma, delta, D) # calc_eigen_value_for_two_site(alpha, beta, gamma, delta, D)
                        if not are_all_real_numbers(eigen_values):
                            print("Not all eigen values are real numbers.")
                            print(f"alpha: {alpha}, beta: {beta}, gamma: {gamma}, delta: {delta}, D: {D}")
                            print(eigen_values)
                            break
                
    print("All eigen values are real numbers.")