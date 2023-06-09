function zi = shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi )

%*****************************************************************************80
%
%% shepard_interp_2d() evaluates a 2D Shepard interpolant.
%
%  Discussion:
%
%    This code should be vectorized.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 August 2012
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Donald Shepard,
%    A two-dimensional interpolation function for irregularly spaced data,
%    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
%    ACM, pages 517-524, 1969.
%
%  Input:
%
%    integer ND, the number of data points.
%
%    real XD(ND,1), YD(ND,1), the data points.
%
%    real ZD(ND,1), the data values.
%
%    real P, the power.
%
%    integer NI, the number of interpolation points.
%
%    real XI(NI,1), YI(NI,1), the interpolation points.
%
%  Output:
%
%    real ZI(NI,1), the interpolated values.
% 
  zi = zeros ( ni, 1 );
%   ni 
%   nd
%   
%   error('er')

  for i = 1 : ni %448 (triangles)

    if ( p == 0.0 )

      w(1:nd,1) = 1.0 / nd;

    else

      w = zeros ( nd, 1 );

      z = -1;
      for j = 1 : nd %930 (collocation points)
        w(j) = sqrt ( ( xi(i) - xd(j) )^2 + ( yi(i) - yd(j) )^2 );
        %[xi(i),xd(j),yi(i),yd(j),w(j)]

        if ( w(j) == 0.0 )
          z = j;
          break
        end
      end

      if ( z ~= -1 )
        w = zeros ( nd, 1 );
        w(z) = 1.0;
      else
        w(1:nd,1) = 1.0 ./ w(1:nd,1) .^ p;
        %================ SUMMATION WITH FOR LOOP
        suma = 0.0;
        for k=1:length(w)
            suma = suma + w(k);
        end
        %=========================================
        s = suma;%sum ( w );
        w(1:nd,1) = w(1:nd,1) / s;
      end

    end
    %================ DOT PRODUCT WITH FOR LOOP
    dotproc = 0.0;
    for k=1:length(w)
        dotproc = dotproc + w(k)*zd(k);
    end
    %=========================================
    %size(w)
    %size(zd)
    zi(i) = dotproc;%w' * zd;

  end

  return
end
