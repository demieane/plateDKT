function zi = shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi )

% https://people.sc.fsu.edu/~jburkardt/c_src/shepard_interp_2d/shepard_interp_2d.html

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

  for i = 1 : ni

    if ( p == 0.0 )

      w(1:nd,1) = 1.0 / nd;

    else

      w = zeros ( nd, 1 );

      z = -1;
      for j = 1 : nd
        w(j) = sqrt ( ( xi(i) - xd(j) )^2 + ( yi(i) - yd(j) )^2 );
        if ( w(j) == 0.0 )
          z = j;
          break
        end
      end

      if ( z ~= -1 )
        w = zeros ( nd, 1 );
        w(z) = 1.0;
      else
        s = 0.0;
        for ii = 1:nd
            w(ii) = 1.0 ./ w(ii,1) .^ p;
            s = s + w(ii);
        end
%         w(1:nd,1) = 1.0 ./ w(1:nd,1) .^ p;
%         s = sum ( w );
        w(1:nd,1) = w(1:nd,1) / s;
      end

    end

%     size(zd)
%     size(w)
    dotproc = 0;
    for k=1:nd
        dotproc = dotproc + w(k)*zd(k);
    end
    zi(i) = dotproc;%w' * zd;

  end
  
  size(zi)

  return
end
