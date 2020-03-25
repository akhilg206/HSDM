%  File: dss004.m
%
   function [ux]=dss004(xl,xu,n,u)
%
%  Function dss004 computes the first derivative, u , of a
%                                                  x
%  variable u over the spatial domain xl le x le xu from classical
%  five-point, fourth-order finite difference approximations
%
%  Argument list
%
%     xl      Lower boundary value of x (input)
%
%     xu      Upper boundary value of x (input)
%
%     n       Number of grid points in the x domain including the
%             boundary points (input)
%
%     u       One-dimensional array containing the values of u at
%             the n grid point points for which the derivative is
%             to be computed (input)
%
%     ux      One-dimensional array containing the numerical
%             values of the derivatives of u at the n grid points
%             (output)
%
%  The mathematical details of the following Taylor series (or
%  polynomials) are given in routine dss002.
%
%  Five-point formulas
%
%     (1)  Left end, point i = 1
%
%                                   2            3            4
%  a(u2 = u1 + u1  ( dx) + u1  ( dx)  + u1  ( dx)  + u1  ( dx)
%                x   1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  ( dx)  + u1  ( dx)  + u1  ( dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                   2            3            4
%  b(u3 = u1 + u1  (2dx) + u1  (2dx)  + u1  (2dx)  + u1  (2dx)
%                x   1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (2dx)  + u1  (2dx)  + u1  (2dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                   2            3            4
%  c(u4 = u1 + u1  (3dx) + u1  (3dx)  + u1  (3dx)  + u1  (3dx)
%                x   1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (3dx)  + u1  (3dx)  + u1  (3dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                   2            3            4
%  d(u5 = u1 + u1  (4dx) + u1  (4dx)  + u1  (4dx)  + u1  (4dx)
%                x   1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u1  (4dx)  + u1  (4dx)  + u1  (4dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%  Constants a, b, c and d are selected so that the coefficients
%  of the u1  terms sum to one and the coefficients of the u1  ,
%           x                                                2x
%  u1   and u1   terms sum to zero
%    3x       4x
%
%  a +   2b +   3c +   4d = 1
%
%  a +   4b +   9c +  16d = 0
%
%  a +   8b +  27c +  64d = 0
%
%  a +  16b +  81c + 256d = 0
%
%  Simultaneous solution for a, b, c and d followed by the solu-
%  tion of the preceding Taylor series, truncated after the u
%                                                            4x
%  terms, for u1  gives the following five-point approximation
%               x
%                                                         4
%  u1  = (1/12dx)(-25u1 + 48u2 - 36u3 + 16u4 - 3u5) + O(dx )   (1)
%    x
%
%     (2)  Interior point, i = 2
%
%                                   2            3            4
%  a(u1 = u2 + u2  (-dx) + u2  (-dx)  + u2  (-dx)  + u2  (-dx)
%                x   1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u2  (-dx)  + u2  (-dx)  + u2  (-dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                   2            3            4
%  b(u3 = u2 + u2  ( dx) + u2  ( dx)  + u2  ( dx)  + u2  ( dx)
%                x   1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u2  ( dx)  + u2  ( dx)  + u2  ( dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                   2            3            4
%  c(u4 = u2 + u2  (2dx) + u2  (2dx)  + u2  (2dx)  + u2  (2dx)
%                x   1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u2  (2dx)  + u2  (2dx)  + u2  (2dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%                                   2            3            4
%  d(u5 = u2 + u2  (3dx) + u2  (3dx)  + u2  (3dx)  + u2  (3dx)
%                x   1f      2x  2f       3x  3f       4x  4f
%
%                      5            6            7
%           + u2  (3dx)  + u2  (3dx)  + u2  (3dx)  + ...)
%               5x  5f       6x  6f       7x  7f
%
%  -a +   b +  2c +  3d = 1
%
%   a +   b +  4c +  9d = 0
%
%  -a +   b +  8c + 27d = 0
%
%   a +   b + 16c + 81d = 0
%
%  Simultaneous solution for a, b, c and d followed by the solu-
%  tion of the preceding Taylor series, truncated after the u
%                                                            4x
%  terms, for u1  gives the following five-point approximation
%               x
%                                                        4
%  u2  = (1/12dx)(-3u1 - 10u2 + 18u3 -  6u4 +  u5) + O(dx )    (2)
%    x
%
%     (3)  Interior point i, i ne 2, n-1
%
%                                        2             3
%  a(ui-2 = ui + ui  (-2dx)  + ui  (-2dx)  + ui  (-2dx)
%                  x    1f       2x   2f       3x   3f
%
%                          4             5             6
%              + ui  (-2dx)  + ui  (-2dx)  + ui  (-2dx)  + ...)
%                  4x   4f       5x   5f       6x   6f
%
%                                        2             3
%  b(ui-1 = ui + ui  ( -dx)  + ui  ( -dx)  + ui  ( -dx)
%                  x    1f       2x   2f       3x   3f
%
%                          4             5             6
%              + ui  ( -dx)  + ui  ( -dx)  + ui  ( -dx)  + ...)
%                  4x   4f       5x   5f       6x   6f
%
%                                        2             3
%  c(ui+1 = ui + ui  (  dx)  + ui  (  dx)  + ui  (  dx)
%                  x    1f       2x   2f       3x   3f
%
%                          4             5             6
%              + ui  (  dx)  + ui  (  dx)  + ui  (  dx)  + ...)
%                  4x   4f       5x   5f       6x   6f
%
%                                        2             3
%  d(ui+2 = ui + ui  ( 2dx)  + ui  ( 2dx)  + ui  ( 2dx)
%                  x    1f       2x   2f       3x   3f
%
%                          4             5             6
%              + ui  ( 2dx)  + ui  ( 2dx)  + ui  ( 2dx)  + ...)
%                  4x   4f       5x   5f       6x   6f
%
%   -2a -   b +   c +  2d = 1
%
%    4a +   b +   c +  4d = 0
%
%   -8a -   b +   c +  8d = 0
%
%   16a +   b +   c + 16d = 0
%
%  Simultaneous solution for a, b, c and d followed by the solu-
%  tion of the preceding Taylor series, truncated after the u
%                                                            4x
%  terms, for u1  gives the following five-point approximation
%               x
%                                                          4
%  ui  = (1/12dx)(ui-2 - 8ui-1 + 0ui + 8ui+1 - ui+2) + O(dx )  (3)
%    x
%
%     (4)  Interior point, i = n-1
%
%                                              2               3
%  a(un-4 = un-1 + un-1  (-3dx)  + un-1  (-3dx)  + un-1  (-3dx)
%                      x    1f         2x   2f         3x   3f
%
%                       4               5               6
%         + un-1  (-3dx)  + un-1  (-3dx)  + un-1  (-3dx)  + ...
%               4x   4f         5x   5f         6x   6f
%
%                                              2               3
%  b(un-3 = un-1 + un-1  (-2dx)  + un-1  (-2dx)  + un-1  (-2dx)
%                      x    1f         2x   2f         3x   3f
%
%                       4               5               6
%         + un-1  (-2dx)  + un-1  (-2dx)  + un-1  (-2dx)  + ...
%               4x   4f         5x   5f         6x   6f
%
%                                              2               3
%  c(un-2 = un-1 + un-1  ( -dx)  + un-1  (- -x)  + un-1  ( -dx)
%                      x    1f         2x   2f         3x   3f
%
%                       4               5               6
%         + un-1  ( -dx)  + un-1  ( -dx)  + un-1  ( -dx)  + ...
%               4x   4f         5x   5f         6x   6f
%
%                                              2               3
%  d(un   = un-1 + un-1  (  dx)  + un-1  (  dx)  + un-1  (  dx)
%                      x    1f         2x   2f         3x   3f
%
%                       4               5               6
%         + un-1  (  dx)  + un-1  (  dx)  + un-1  (  dx)  + ...
%               4x   4f         5x   5f         6x   6f
%
%  -3a -  2b -   c +   d = 1
%
%   9a +  4b +   c +   d = 0
%
% -27a -  8b -   c +   d = 0
%
%  81a + 16b +   c +   d = 0
%
%  Simultaneous solution for a, b, c and d followed by the solu-
%  tion of the preceding Taylor series, truncated after the u
%                                                            4x
%  terms, for u1  gives the following five-point approximation
%               x
%                                                                4
%  un-1  = (1/12dx)(-un-4 + 6un-3 - 18un-2 + 10un-1 + 3un) + O(dx )
%      x
%                                                              (4)
%
%    (5)  Right end, point i = n
%
%                                       2             3
%  a(un-4 = un + un (-4dx)  + un  (-4dx)  + un  (-4dx)
%                  x   1f       2x   2f       3x   3f
%
%                         4             5             6
%             + un  (-4dx)  + un  (-4dx)  + un  (-4dx)  + ...)
%                 4x   4f       5x   5f       6x   6f
%
%                                       2             3
%  b(un-3 = un + un (-3dx)  + un  (-3dx)  + un  (-3dx)
%                  x   1f       2x   2f       3x   3f
%
%                         4             5             6
%             + un  (-3dx)  + un  (-3dx)  + un  (-3dx)  + ...)
%                 4x   4f       5x   5f       6x   6f
%
%                                       2             3
%  c(un-2 = un + un (-2dx)  + un  (-2dx)  + un  (-2dx)
%                  x   1f       2x   2f       3x   3f
%
%                         4             5             6
%             + un  (-2dx)  + un  (-2dx)  + un  (-2dx)  + ...)
%                 4x   4f       5x   5f       6x   6f
%
%                                       2             3
%  d(un-1 = un + un ( -dx)  + un  ( -dx)  + un  ( -dx)
%                  x   1f       2x   2f       3x   3f
%
%                         4             5             6
%             + un  ( -dx)  + un  ( -dx)  + un  ( -dx)  + ...)
%                 4x   4f       5x   5f       6x   6f
%
%   -4a -  3b -  2c -   d = 1
%
%   16a +  9b +  4c +   d = 0
%
%  -64a - 27b -  8c -   d = 0
%
%  256a + 81b + 16c +   d = 0
%
%  Simultaneous solution for a, b, c and d followed by the solu-
%  tion of the preceding Taylor series, truncated after the u
%                                                            4x
%  terms, for u1  gives the following five-point approximation
%               x
%                                                                4
%  un  = (1/12dx)(3un-4 - 16un-3 + 36un-2 - 48un-1 + 25un) + O(dx )
%    x
%                                                              (5)
%
%  The weighting coefficients for equations (1) to (5) can be
%  summarized as
%
%             -25   48  -36   16   -3
%
%              -3  -10   18   -6    1
%
%       1/12    1   -8    0    8   -1
%
%              -1    6  -18   10    3
%
%               3  -16   36  -48   25
%
%  which are the coefficients reported by Bickley for n = 4, m =
%  1, p = 0, 1, 2, 3, 4 (Bickley, W. G., Formulae for Numerical
%  Differentiation, Math. Gaz., vol. 25, 1941.  Note - the Bickley
%  coefficients have been divided by a common factor of two).
%
%  Equations (1) to (5) can now be programmed to generate the
%  derivative u (x) of function u(x).
%              x
%
%  Compute the spatial increment
   dx=(xu-xl)/(n-1);
   r4fdx=1./(12.*dx);
   nm2=n-2;
%
%  Equation (1) (note - the rhs of equations (1), (2), (3), (4)
%  and (5) have been formatted so that the numerical weighting
%  coefficients can be more easily associated with the Bickley
%  matrix above)
   ux(  1)=r4fdx*...
     ( -25.*u(  1) +48.*u(  2) -36.*u(  3) +16.*u(  4)  -3.*u(  5));
%
%  Equation (2)
   ux(  2)=r4fdx*...
     (  -3.*u(  1) -10.*u(  2) +18.*u(  3)  -6.*u(  4)  +1.*u(  5));
%
%  Equation (3)
   for i=3:nm2
     ux(  i)=r4fdx*...
     (  +1.*u(i-2)  -8.*u(i-1)  +0.*u(  i)  +8.*u(i+1)  -1.*u(i+2));
   end
%
%  Equation (4)
   ux(n-1)=r4fdx*...
     (  -1.*u(n-4)  +6.*u(n-3) -18.*u(n-2) +10.*u(n-1)  +3.*u(  n));
%
%  Equation (5)
   ux(  n)=r4fdx*...
     (   3.*u(n-4) -16.*u(n-3) +36.*u(n-2) -48.*u(n-1) +25.*u(  n));
