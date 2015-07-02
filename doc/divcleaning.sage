sage: nx,ny,nz = var('nx, ny, nz')
sage: (nx > 0).assume()
sage: (ny > 0).assume()
sage: (nz > 0).assume()
sage: n = matrix([nx, ny, nz])
sage: nt = n.transpose()
sage: c1, c2 = var('c1, c2')
sage: (c1 > 0).assume()
sage: (c2 > 0).assume()
sage: rot = matrix([[0,-nz,ny],[nz,0,-nx],[-ny,nx,0]])
sage: r = var('r')
sage: (r > 0).assume()
sage: An = block_matrix([[0*rot,-rot,c1*nt,0*nt],[rot,0*rot,0*nt,c2*nt],[c1*n,0*n,0,0],[0*n,c2*n,0,0]], subdivide=False)
sage: An
[    0     0     0     0    nz   -ny c1*nx     0]
[    0     0     0   -nz     0    nx c1*ny     0]
[    0     0     0    ny   -nx     0 c1*nz     0]
[    0   -nz    ny     0     0     0     0 c2*nx]
[   nz     0   -nx     0     0     0     0 c2*ny]
[  -ny    nx     0     0     0     0     0 c2*nz]
[c1*nx c1*ny c1*nz     0     0     0     0     0]
[    0     0     0 c2*nx c2*ny c2*nz     0     0]
sage: eigv = An.eigenvalues()
sage: eigv
[-sqrt(nx^2 + ny^2 + nz^2),
 -sqrt(nx^2 + ny^2 + nz^2),
 sqrt(nx^2 + ny^2 + nz^2),
 sqrt(nx^2 + ny^2 + nz^2),
 -sqrt(nx^2 + ny^2 + nz^2)*c1,
 sqrt(nx^2 + ny^2 + nz^2)*c1,
 -sqrt(nx^2 + ny^2 + nz^2)*c2,
 sqrt(nx^2 + ny^2 + nz^2)*c2]
sage: [v.subs({sqrt(nx*nx + ny * ny + nz * nz) : r}) for v in eigv]
[-r, -r, r, r, -c1*r, c1*r, -c2*r, c2*r]
sage: eigV = An.eigenvectors_right()
sage: eigV
[(-sqrt(nx^2 + ny^2 + nz^2),
  [(1, 0, -nx/nz, sqrt(nx^2 + ny^2 + nz^2)*nx*ny/(nz^3 + (nx^2 + ny^2)*nz), -sqrt(nx^2 + ny^2 + nz^2)*(nx^2 + nz^2)/(nz^3 + (nx^2 + ny^2)*nz), ny/sqrt(nx^2 + ny^2 + nz^2), 0, 0),
   (0, 1, -ny/nz, sqrt(nx^2 + ny^2 + nz^2)*(ny^2 + nz^2)/(nz^3 + (nx^2 + ny^2)*nz), -sqrt(nx^2 + ny^2 + nz^2)*nx*ny/(nz^3 + (nx^2 + ny^2)*nz), -nx/sqrt(nx^2 + ny^2 + nz^2), 0, 0)],
  2),
 (sqrt(nx^2 + ny^2 + nz^2),
  [(1, 0, -nx/nz, -sqrt(nx^2 + ny^2 + nz^2)*nx*ny/(nz^3 + (nx^2 + ny^2)*nz), sqrt(nx^2 + ny^2 + nz^2)*(nx^2 + nz^2)/(nz^3 + (nx^2 + ny^2)*nz), -ny/sqrt(nx^2 + ny^2 + nz^2), 0, 0),
   (0, 1, -ny/nz, -sqrt(nx^2 + ny^2 + nz^2)*(ny^2 + nz^2)/(nz^3 + (nx^2 + ny^2)*nz), sqrt(nx^2 + ny^2 + nz^2)*nx*ny/(nz^3 + (nx^2 + ny^2)*nz), nx/sqrt(nx^2 + ny^2 + nz^2), 0, 0)],
  2),
 (-sqrt(nx^2 + ny^2 + nz^2)*c1,
  [(1, ny/nx, nz/nx, 0, 0, 0, -sqrt(nx^2 + ny^2 + nz^2)/nx, 0)],
  1),
 (sqrt(nx^2 + ny^2 + nz^2)*c1,
  [(1, ny/nx, nz/nx, 0, 0, 0, sqrt(nx^2 + ny^2 + nz^2)/nx, 0)],
  1),
 (-sqrt(nx^2 + ny^2 + nz^2)*c2,
  [(0, 0, 0, 1, ny/nx, nz/nx, 0, -sqrt(nx^2 + ny^2 + nz^2)/nx)],
  1),
 (sqrt(nx^2 + ny^2 + nz^2)*c2,
  [(0, 0, 0, 1, ny/nx, nz/nx, 0, sqrt(nx^2 + ny^2 + nz^2)/nx)],
  1)]
sage: eigV[0][1][0]
(1, 0, -nx/nz, sqrt(nx^2 + ny^2 + nz^2)*nx*ny/(nz^3 + (nx^2 + ny^2)*nz), -sqrt(nx^2 + ny^2 + nz^2)*(nx^2 + nz^2)/(nz^3 + (nx^2 + ny^2)*nz), ny/sqrt(nx^2 + ny^2 + nz^2), 0, 0)
sage: D, P = An.eigenmatrix_right()
sage: factor(expand(An*P)) == factor(expand(P*D))
True
sage: Pinv = factor(expand(P^(-1)))
sage: Pinv
[1/2*(ny^2 + nz^2)/(nx^2 + ny^2 + nz^2)        -1/2*nx*ny/(nx^2 + ny^2 + nz^2)        -1/2*nx*nz/(nx^2 + ny^2 + nz^2)                                      0       -1/2*nz/sqrt(nx^2 + ny^2 + nz^2)        1/2*ny/sqrt(nx^2 + ny^2 + nz^2)                                      0                                      0]
[       -1/2*nx*ny/(nx^2 + ny^2 + nz^2) 1/2*(nx^2 + nz^2)/(nx^2 + ny^2 + nz^2)        -1/2*ny*nz/(nx^2 + ny^2 + nz^2)        1/2*nz/sqrt(nx^2 + ny^2 + nz^2)                                      0       -1/2*nx/sqrt(nx^2 + ny^2 + nz^2)                                      0                                      0]
[1/2*(ny^2 + nz^2)/(nx^2 + ny^2 + nz^2)        -1/2*nx*ny/(nx^2 + ny^2 + nz^2)        -1/2*nx*nz/(nx^2 + ny^2 + nz^2)                                      0        1/2*nz/sqrt(nx^2 + ny^2 + nz^2)       -1/2*ny/sqrt(nx^2 + ny^2 + nz^2)                                      0                                      0]
[       -1/2*nx*ny/(nx^2 + ny^2 + nz^2) 1/2*(nx^2 + nz^2)/(nx^2 + ny^2 + nz^2)        -1/2*ny*nz/(nx^2 + ny^2 + nz^2)       -1/2*nz/sqrt(nx^2 + ny^2 + nz^2)                                      0        1/2*nx/sqrt(nx^2 + ny^2 + nz^2)                                      0                                      0]
[         1/2*nx^2/(nx^2 + ny^2 + nz^2)         1/2*nx*ny/(nx^2 + ny^2 + nz^2)         1/2*nx*nz/(nx^2 + ny^2 + nz^2)                                      0                                      0                                      0       -1/2*nx/sqrt(nx^2 + ny^2 + nz^2)                                      0]
[         1/2*nx^2/(nx^2 + ny^2 + nz^2)         1/2*nx*ny/(nx^2 + ny^2 + nz^2)         1/2*nx*nz/(nx^2 + ny^2 + nz^2)                                      0                                      0                                      0        1/2*nx/sqrt(nx^2 + ny^2 + nz^2)                                      0]
[                                     0                                      0                                      0          1/2*nx^2/(nx^2 + ny^2 + nz^2)         1/2*nx*ny/(nx^2 + ny^2 + nz^2)         1/2*nx*nz/(nx^2 + ny^2 + nz^2)                                      0       -1/2*nx/sqrt(nx^2 + ny^2 + nz^2)]
[                                     0                                      0                                      0          1/2*nx^2/(nx^2 + ny^2 + nz^2)         1/2*nx*ny/(nx^2 + ny^2 + nz^2)         1/2*nx*nz/(nx^2 + ny^2 + nz^2)                                      0        1/2*nx/sqrt(nx^2 + ny^2 + nz^2)]
sage: Pinv.subs({nx^2 + ny^2 + nz^2 : r^2})
[1/2*(ny^2 + nz^2)/r^2        -1/2*nx*ny/r^2        -1/2*nx*nz/r^2                     0     -1/2*nz/sqrt(r^2)      1/2*ny/sqrt(r^2)                     0                     0]
[       -1/2*nx*ny/r^2 1/2*(nx^2 + nz^2)/r^2        -1/2*ny*nz/r^2      1/2*nz/sqrt(r^2)                     0     -1/2*nx/sqrt(r^2)                     0                     0]
[1/2*(ny^2 + nz^2)/r^2        -1/2*nx*ny/r^2        -1/2*nx*nz/r^2                     0      1/2*nz/sqrt(r^2)     -1/2*ny/sqrt(r^2)                     0                     0]
[       -1/2*nx*ny/r^2 1/2*(nx^2 + nz^2)/r^2        -1/2*ny*nz/r^2     -1/2*nz/sqrt(r^2)                     0      1/2*nx/sqrt(r^2)                     0                     0]
[         1/2*nx^2/r^2         1/2*nx*ny/r^2         1/2*nx*nz/r^2                     0                     0                     0     -1/2*nx/sqrt(r^2)                     0]
[         1/2*nx^2/r^2         1/2*nx*ny/r^2         1/2*nx*nz/r^2                     0                     0                     0      1/2*nx/sqrt(r^2)                     0]
[                    0                     0                     0          1/2*nx^2/r^2         1/2*nx*ny/r^2         1/2*nx*nz/r^2                     0     -1/2*nx/sqrt(r^2)]
[                    0                     0                     0          1/2*nx^2/r^2         1/2*nx*ny/r^2         1/2*nx*nz/r^2                     0      1/2*nx/sqrt(r^2)]
sage: D.subs({sqrt(nx*nx + ny * ny + nz * nz) : r})
[   -r     0     0     0     0     0     0     0]
[    0    -r     0     0     0     0     0     0]
[    0     0     r     0     0     0     0     0]
[    0     0     0     r     0     0     0     0]
[    0     0     0     0 -c1*r     0     0     0]
[    0     0     0     0     0  c1*r     0     0]
[    0     0     0     0     0     0 -c2*r     0]
[    0     0     0     0     0     0     0  c2*r]
sage: pD = 1 * D
sage: i = 0
sage: while i < 8:
...       if pD[i,i] < 0:
...           pD[i,i] = 0
...       i += 1
sage: nD = 1 * D
sage: i = 0
sage: while i < 8:
...       if nD[i,i] > 0:
...           nD[i,i] = 0
...       i += 1
sage: negA = factor(expand(P * nD * Pinv))
sage: nicenegA = negA.subs({1/sqrt(nx^2 + ny^2 + nz^2): 1/r})
sage: nicenegA
[  -1/2*(c1*nx^2 + ny^2 + nz^2)/r            -1/2*(c1 - 1)*nx*ny/r            -1/2*(c1 - 1)*nx*nz/r                                0                           1/2*nz                          -1/2*ny                        1/2*c1*nx                                0]
[           -1/2*(c1 - 1)*nx*ny/r   -1/2*(c1*ny^2 + nx^2 + nz^2)/r            -1/2*(c1 - 1)*ny*nz/r                          -1/2*nz                                0                           1/2*nx                        1/2*c1*ny                                0]
[           -1/2*(c1 - 1)*nx*nz/r            -1/2*(c1 - 1)*ny*nz/r   -1/2*(c1*nz^2 + nx^2 + ny^2)/r                           1/2*ny                          -1/2*nx                                0                        1/2*c1*nz                                0]
[                               0                          -1/2*nz                           1/2*ny   -1/2*(c2*nx^2 + ny^2 + nz^2)/r            -1/2*(c2 - 1)*nx*ny/r            -1/2*(c2 - 1)*nx*nz/r                                0                        1/2*c2*nx]
[                          1/2*nz                                0                          -1/2*nx            -1/2*(c2 - 1)*nx*ny/r   -1/2*(c2*ny^2 + nx^2 + nz^2)/r            -1/2*(c2 - 1)*ny*nz/r                                0                        1/2*c2*ny]
[                         -1/2*ny                           1/2*nx                                0            -1/2*(c2 - 1)*nx*nz/r            -1/2*(c2 - 1)*ny*nz/r   -1/2*(c2*nz^2 + nx^2 + ny^2)/r                                0                        1/2*c2*nz]
[                       1/2*c1*nx                        1/2*c1*ny                        1/2*c1*nz                                0                                0                                0 -1/2*sqrt(nx^2 + ny^2 + nz^2)*c1                                0]
[                               0                                0                                0                        1/2*c2*nx                        1/2*c2*ny                        1/2*c2*nz                                0 -1/2*sqrt(nx^2 + ny^2 + nz^2)*c2]
sage: latex(nicenegA)
\left(\begin{array}{rrrrrrrr}
-\frac{c_{1} \mathit{nx}^{2} + \mathit{ny}^{2} + \mathit{nz}^{2}}{2 \, r} & -\frac{{\left(c_{1} - 1\right)} \mathit{nx} \mathit{ny}}{2 \, r} & -\frac{{\left(c_{1} - 1\right)} \mathit{nx} \mathit{nz}}{2 \, r} & 0 & \frac{1}{2} \, \mathit{nz} & -\frac{1}{2} \, \mathit{ny} & \frac{1}{2} \, c_{1} \mathit{nx} & 0 \\
-\frac{{\left(c_{1} - 1\right)} \mathit{nx} \mathit{ny}}{2 \, r} & -\frac{c_{1} \mathit{ny}^{2} + \mathit{nx}^{2} + \mathit{nz}^{2}}{2 \, r} & -\frac{{\left(c_{1} - 1\right)} \mathit{ny} \mathit{nz}}{2 \, r} & -\frac{1}{2} \, \mathit{nz} & 0 & \frac{1}{2} \, \mathit{nx} & \frac{1}{2} \, c_{1} \mathit{ny} & 0 \\
-\frac{{\left(c_{1} - 1\right)} \mathit{nx} \mathit{nz}}{2 \, r} & -\frac{{\left(c_{1} - 1\right)} \mathit{ny} \mathit{nz}}{2 \, r} & -\frac{c_{1} \mathit{nz}^{2} + \mathit{nx}^{2} + \mathit{ny}^{2}}{2 \, r} & \frac{1}{2} \, \mathit{ny} & -\frac{1}{2} \, \mathit{nx} & 0 & \frac{1}{2} \, c_{1} \mathit{nz} & 0 \\
0 & -\frac{1}{2} \, \mathit{nz} & \frac{1}{2} \, \mathit{ny} & -\frac{c_{2} \mathit{nx}^{2} + \mathit{ny}^{2} + \mathit{nz}^{2}}{2 \, r} & -\frac{{\left(c_{2} - 1\right)} \mathit{nx} \mathit{ny}}{2 \, r} & -\frac{{\left(c_{2} - 1\right)} \mathit{nx} \mathit{nz}}{2 \, r} & 0 & \frac{1}{2} \, c_{2} \mathit{nx} \\
\frac{1}{2} \, \mathit{nz} & 0 & -\frac{1}{2} \, \mathit{nx} & -\frac{{\left(c_{2} - 1\right)} \mathit{nx} \mathit{ny}}{2 \, r} & -\frac{c_{2} \mathit{ny}^{2} + \mathit{nx}^{2} + \mathit{nz}^{2}}{2 \, r} & -\frac{{\left(c_{2} - 1\right)} \mathit{ny} \mathit{nz}}{2 \, r} & 0 & \frac{1}{2} \, c_{2} \mathit{ny} \\
-\frac{1}{2} \, \mathit{ny} & \frac{1}{2} \, \mathit{nx} & 0 & -\frac{{\left(c_{2} - 1\right)} \mathit{nx} \mathit{nz}}{2 \, r} & -\frac{{\left(c_{2} - 1\right)} \mathit{ny} \mathit{nz}}{2 \, r} & -\frac{c_{2} \mathit{nz}^{2} + \mathit{nx}^{2} + \mathit{ny}^{2}}{2 \, r} & 0 & \frac{1}{2} \, c_{2} \mathit{nz} \\
\frac{1}{2} \, c_{1} \mathit{nx} & \frac{1}{2} \, c_{1} \mathit{ny} & \frac{1}{2} \, c_{1} \mathit{nz} & 0 & 0 & 0 & -\frac{1}{2} \, \sqrt{\mathit{nx}^{2} + \mathit{ny}^{2} + \mathit{nz}^{2}} c_{1} & 0 \\
0 & 0 & 0 & \frac{1}{2} \, c_{2} \mathit{nx} & \frac{1}{2} \, c_{2} \mathit{ny} & \frac{1}{2} \, c_{2} \mathit{nz} & 0 & -\frac{1}{2} \, \sqrt{\mathit{nx}^{2} + \mathit{ny}^{2} + \mathit{nz}^{2}} c_{2}
\end{array}\right)
sage: posA = factor(expand(P * pD * Pinv))
sage: niceposA = posA.subs({1/sqrt(nx^2 + ny^2 + nz^2): 1/r})
sage: niceposA
[  1/2*(c1*nx^2 + ny^2 + nz^2)/r            1/2*(c1 - 1)*nx*ny/r            1/2*(c1 - 1)*nx*nz/r                               0                          1/2*nz                         -1/2*ny                       1/2*c1*nx                               0]
[           1/2*(c1 - 1)*nx*ny/r   1/2*(c1*ny^2 + nx^2 + nz^2)/r            1/2*(c1 - 1)*ny*nz/r                         -1/2*nz                               0                          1/2*nx                       1/2*c1*ny                               0]
[           1/2*(c1 - 1)*nx*nz/r            1/2*(c1 - 1)*ny*nz/r   1/2*(c1*nz^2 + nx^2 + ny^2)/r                          1/2*ny                         -1/2*nx                               0                       1/2*c1*nz                               0]
[                              0                         -1/2*nz                          1/2*ny   1/2*(c2*nx^2 + ny^2 + nz^2)/r            1/2*(c2 - 1)*nx*ny/r            1/2*(c2 - 1)*nx*nz/r                               0                       1/2*c2*nx]
[                         1/2*nz                               0                         -1/2*nx            1/2*(c2 - 1)*nx*ny/r   1/2*(c2*ny^2 + nx^2 + nz^2)/r            1/2*(c2 - 1)*ny*nz/r                               0                       1/2*c2*ny]
[                        -1/2*ny                          1/2*nx                               0            1/2*(c2 - 1)*nx*nz/r            1/2*(c2 - 1)*ny*nz/r   1/2*(c2*nz^2 + nx^2 + ny^2)/r                               0                       1/2*c2*nz]
[                      1/2*c1*nx                       1/2*c1*ny                       1/2*c1*nz                               0                               0                               0 1/2*sqrt(nx^2 + ny^2 + nz^2)*c1                               0]
[                              0                               0                               0                       1/2*c2*nx                       1/2*c2*ny                       1/2*c2*nz                               0 1/2*sqrt(nx^2 + ny^2 + nz^2)*c2]
sage: latex(niceposA)
\left(\begin{array}{rrrrrrrr}
\frac{c_{1} \mathit{nx}^{2} + \mathit{ny}^{2} + \mathit{nz}^{2}}{2 \, r} & \frac{{\left(c_{1} - 1\right)} \mathit{nx} \mathit{ny}}{2 \, r} & \frac{{\left(c_{1} - 1\right)} \mathit{nx} \mathit{nz}}{2 \, r} & 0 & \frac{1}{2} \, \mathit{nz} & -\frac{1}{2} \, \mathit{ny} & \frac{1}{2} \, c_{1} \mathit{nx} & 0 \\
\frac{{\left(c_{1} - 1\right)} \mathit{nx} \mathit{ny}}{2 \, r} & \frac{c_{1} \mathit{ny}^{2} + \mathit{nx}^{2} + \mathit{nz}^{2}}{2 \, r} & \frac{{\left(c_{1} - 1\right)} \mathit{ny} \mathit{nz}}{2 \, r} & -\frac{1}{2} \, \mathit{nz} & 0 & \frac{1}{2} \, \mathit{nx} & \frac{1}{2} \, c_{1} \mathit{ny} & 0 \\
\frac{{\left(c_{1} - 1\right)} \mathit{nx} \mathit{nz}}{2 \, r} & \frac{{\left(c_{1} - 1\right)} \mathit{ny} \mathit{nz}}{2 \, r} & \frac{c_{1} \mathit{nz}^{2} + \mathit{nx}^{2} + \mathit{ny}^{2}}{2 \, r} & \frac{1}{2} \, \mathit{ny} & -\frac{1}{2} \, \mathit{nx} & 0 & \frac{1}{2} \, c_{1} \mathit{nz} & 0 \\
0 & -\frac{1}{2} \, \mathit{nz} & \frac{1}{2} \, \mathit{ny} & \frac{c_{2} \mathit{nx}^{2} + \mathit{ny}^{2} + \mathit{nz}^{2}}{2 \, r} & \frac{{\left(c_{2} - 1\right)} \mathit{nx} \mathit{ny}}{2 \, r} & \frac{{\left(c_{2} - 1\right)} \mathit{nx} \mathit{nz}}{2 \, r} & 0 & \frac{1}{2} \, c_{2} \mathit{nx} \\
\frac{1}{2} \, \mathit{nz} & 0 & -\frac{1}{2} \, \mathit{nx} & \frac{{\left(c_{2} - 1\right)} \mathit{nx} \mathit{ny}}{2 \, r} & \frac{c_{2} \mathit{ny}^{2} + \mathit{nx}^{2} + \mathit{nz}^{2}}{2 \, r} & \frac{{\left(c_{2} - 1\right)} \mathit{ny} \mathit{nz}}{2 \, r} & 0 & \frac{1}{2} \, c_{2} \mathit{ny} \\
-\frac{1}{2} \, \mathit{ny} & \frac{1}{2} \, \mathit{nx} & 0 & \frac{{\left(c_{2} - 1\right)} \mathit{nx} \mathit{nz}}{2 \, r} & \frac{{\left(c_{2} - 1\right)} \mathit{ny} \mathit{nz}}{2 \, r} & \frac{c_{2} \mathit{nz}^{2} + \mathit{nx}^{2} + \mathit{ny}^{2}}{2 \, r} & 0 & \frac{1}{2} \, c_{2} \mathit{nz} \\
\frac{1}{2} \, c_{1} \mathit{nx} & \frac{1}{2} \, c_{1} \mathit{ny} & \frac{1}{2} \, c_{1} \mathit{nz} & 0 & 0 & 0 & \frac{1}{2} \, \sqrt{\mathit{nx}^{2} + \mathit{ny}^{2} + \mathit{nz}^{2}} c_{1} & 0 \\
0 & 0 & 0 & \frac{1}{2} \, c_{2} \mathit{nx} & \frac{1}{2} \, c_{2} \mathit{ny} & \frac{1}{2} \, c_{2} \mathit{nz} & 0 & \frac{1}{2} \, \sqrt{\mathit{nx}^{2} + \mathit{ny}^{2} + \mathit{nz}^{2}} c_{2}
\end{array}\right)
sage: An
[    0     0     0     0    nz   -ny c1*nx     0]
[    0     0     0   -nz     0    nx c1*ny     0]
[    0     0     0    ny   -nx     0 c1*nz     0]
[    0   -nz    ny     0     0     0     0 c2*nx]
[   nz     0   -nx     0     0     0     0 c2*ny]
[  -ny    nx     0     0     0     0     0 c2*nz]
[c1*nx c1*ny c1*nz     0     0     0     0     0]
[    0     0     0 c2*nx c2*ny c2*nz     0     0]
sage: latex(An)
\left(\begin{array}{rrrrrrrr}
0 & 0 & 0 & 0 & \mathit{nz} & -\mathit{ny} & c_{1} \mathit{nx} & 0 \\
0 & 0 & 0 & -\mathit{nz} & 0 & \mathit{nx} & c_{1} \mathit{ny} & 0 \\
0 & 0 & 0 & \mathit{ny} & -\mathit{nx} & 0 & c_{1} \mathit{nz} & 0 \\
0 & -\mathit{nz} & \mathit{ny} & 0 & 0 & 0 & 0 & c_{2} \mathit{nx} \\
\mathit{nz} & 0 & -\mathit{nx} & 0 & 0 & 0 & 0 & c_{2} \mathit{ny} \\
-\mathit{ny} & \mathit{nx} & 0 & 0 & 0 & 0 & 0 & c_{2} \mathit{nz} \\
c_{1} \mathit{nx} & c_{1} \mathit{ny} & c_{1} \mathit{nz} & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & c_{2} \mathit{nx} & c_{2} \mathit{ny} & c_{2} \mathit{nz} & 0 & 0
\end{array}\right)
