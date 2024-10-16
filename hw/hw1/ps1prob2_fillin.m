%--------------------------------------------------------------------------
%Matlab script to calculate derivative using three different finite 
%difference formulas
%--------------------------------------------------------------------------
%Problem Set #1, Problem 3
%--------------------------------------------------------------------------
%Author: Guorui Wei (危国锐) (313017602@qq.com)
%Date created: Oct. 17, 2024
%--------------------------------------------------------------------------

close all %close all figure windows
clear all %clear memory of all variables

%-------------------------------------------------------
%Define variables that do not change for different grids
%-------------------------------------------------------
L=3; %Length of domain
m=1; %counter for number of grid sizes to plot

%---------------------------------------------------------------
%Loop over different grid sizes (N = number of points in domain)
%---------------------------------------------------------------
for N = [8 16 32 64 128 256 512 1024 2048] 

  %-------------------------------------
  %Define resolution-dependent variables
  %-------------------------------------
  deltax = L/N; %Grid spacing
  x =0:deltax:L; %define x, grid
  nx=length(x); %number of points in grid

  f = sin(5*x); %function to be differentiated

  dfdx_exact = 5 * cos(5*x); %exact (analytical) derivative

  %-----------------------------------------------------------------------
  %Initialize arrays for 3 schemes with NaNs so that un-used values at the 
  %edge of the domain will not be plotted. This isn't necessary, but just
  %makes the plotting a bit easier.
  %-----------------------------------------------------------------------
  dfdx_1st = NaN*ones(size(dfdx_exact));
  dfdx_2nd = NaN*ones(size(dfdx_exact));
  dfdx_4th = NaN*ones(size(dfdx_exact));

  %-----------------------------------------------------------------------
  %Compute derivatives with three different approximations
  %-----------------------------------------------------------------------
  for i=3:nx-2 %Notice the limits here
      dfdx_1st(i) =  (f(i+1) - f(i)) / deltax; %first order forward
      dfdx_2nd(i) =  (f(i+1) - f(i-1)) / 2 / deltax; %second order central
      dfdx_4th(i) =  (f(i-2) - 8*f(i-1) + 8*f(i+1) - f(i+2)) / 12 / deltax; %fourth order central
  end

  %-----------------------------------------------------------------------
  %Plot the exact and approximated derivatives for the case where N=16
  %-----------------------------------------------------------------------
  if (N == 16)
    figure
    plot(x,dfdx_exact,'-',x,dfdx_1st,'-.',x,dfdx_2nd,'--',x,dfdx_4th,':')
    xlabel('\it x')
    ylabel('d{\itf} / d{\itx}')
    legend('Exact','Forward difference - 1st order',...
          'Central difference - 2nd order', 'Central difference - 4th order')
    set(gca, FontName="Times New Roman", FontSize=10.5)
    print(gcf, "fig_3_1_first_diff_compare.svg", "-dsvg")
  end

  %-----------------------------------------------------------------------
  %Extract the errors for x = 1.5, the midpoint of the domain
  %-----------------------------------------------------------------------
  dx(m) = deltax; %store the current deltax in an array for plotting later
  midpoint = floor((1+nx)/2); %calculate the midpoint index 
  error_1st(m) = abs(dfdx_1st(midpoint)-dfdx_exact(midpoint));
  error_2nd(m) = abs(dfdx_2nd(midpoint)-dfdx_exact(midpoint));
  error_4th(m) = abs(dfdx_4th(midpoint)-dfdx_exact(midpoint));

  m = m +1; %Advance counter for number of grids
  
end %end of loop over different grid sizes

%-----------------------------
%Calculate slopes of lines
%-----------------------------
slope_1st = diff(log(error_1st))./diff(log(dx))
slope_2nd = diff(log(error_2nd))./diff(log(dx))
slope_4th = diff(log(error_4th))./diff(log(dx))

%--------------------------------
%Plot dx versus error at x = 1.5
%--------------------------------
figure
loglog(dx,error_1st,'-*',dx,error_2nd,'x--',dx,error_4th,'+:')
xlabel('\Delta')
ylabel('Error at{\it x} = 1.5')
legend('Forward difference - 1st order','Central difference - 2nd order', ...
      'Central difference - 4th order', Location='southeast')
set(gca, FontName="Times New Roman", FontSize=10.5)
print(gcf, "fig_3_2_first_diff_error_loglog.svg", "-dsvg")
