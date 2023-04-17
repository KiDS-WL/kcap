import mpmath

def roots_normalization(n_max, z_max, return_float=True, 
                        dps=40, J_dps=130, c_dps=80, N_dps=50):
    """Calculate the roots and normalizations for the log-COSEBIS."""
    n_max_idx = n_max+1
    
    r = {}
    N = {}

    with mpmath.workdps(dps):
        # Calculate J's (eq. 32)
        with mpmath.workdps(J_dps):
            J = {k : mpmath.matrix([-mpmath.gammainc(j+1, b=0, a=-k*z_max)/(-k)**(j+1) for j in range(2*n_max_idx+1)]) for k in [1, 2, 4]}

        # Calculate c_10, c_11, c_12
        c = {1 : [(-J[2][2] + J[4][2]/J[4][1]*J[2][1]) / (J[2][0] - J[4][0]/J[4][1]*J[2][1]), #c_10
                  (-J[2][2] + J[4][2]/J[4][0]*J[2][0]) / (J[2][1] - J[4][1]/J[4][0]*J[2][0]), #c_11
                  1] #c_12 
            }

        for n in range(2, n_max_idx):
            n_idx = n-1

            A = mpmath.matrix(n+1,n+1)
            B = mpmath.ones(n+1,1)
            B[n_idx] = -1
            B[n_idx+1] = -1

            # Calculate A
            for j in range(0, n+1): # j=0,...,n
                A[n_idx,j] = J[2][j]/J[2][n+1]   # A_nj
                A[n_idx+1,j] = J[4][j]/J[4][n+1] # A_(n+1)j

                for m in range(1, n): # m=1,...,n-1
                    m_idx = m-1
                    # A_mj
                    A[m_idx,j] = - sum([J[1][i+j]*c[m][i] for i in range(0,m+2)]) \
                                    / sum([J[1][i+n+1]*c[m][i] for i in range(0,m+2)]) #i=0,...,m+1

            with mpmath.workdps(c_dps):
                C = mpmath.lu_solve(A, B)
                c[n] = list(C) + [1] #c_n(n+1) = 1
        
        # Calculate roots and normalization
        for n in c.keys():
            r[n], err = mpmath.polyroots(c[n][::-1], error=True, maxsteps=100, extraprec=30)

            with mpmath.workdps(N_dps):
                t = lambda z: mpmath.polyval(c[n][::-1], z)**2 * mpmath.exp(z)
                N[n] = mpmath.sqrt((mpmath.exp(z_max)-1)/mpmath.quad(t, (0, z_max)))
                
            if return_float:
                r[n] = [float(_) for _ in r[n]]
                N[n] = float(N[n])
                
    return r, N