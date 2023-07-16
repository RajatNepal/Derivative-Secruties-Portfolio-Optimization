% Nelson-Siegel Yield Curve Fitting
% Example with realistic values

% Define the Nelson-Siegel Function
NS = @(b0,b1,b2,tau,t) b0 + b1*(1-exp(-t/tau)).*(t/tau) + b2*(1-exp(-t/tau)).*(t/tau)-b2*exp(-t/tau);

% Define bond data
bond_maturities = [1 2 3 5 7 10 20 30];
bond_yield = [0.023 0.034 0.045 0.058 0.065 0.072 0.087 0.093];

% Fit the Nelson-Siegel model to the bond data
% Starting values for parameters
start = [0.05, -0.02, 0.01, 1];
% Optimization function
f = @(params) sum((NS(params(1),params(2),params(3),params(4),bond_maturities)-bond_yield).^2);
% Optimization
params = fminsearch(f, start);

% Calculate the Yield Curve using the Nelson-Siegel Model
NS_Yield_Curve = NS(params(1),params(2),params(3),params(4),bond_maturities);

% Plot the Yield Curve
plot(bond_maturities, bond_yield, 'ro', bond_maturities, NS_Yield_Curve, 'b-')
xlabel('Maturity')
ylabel('Yield')
legend('Market Yield', 'Nelson-Siegel Yield Curve')


% Select yield for a specific time to maturity
time_to_maturity = 98/12;
yield_at_maturity = interp1(bond_maturities, NS_Yield_Curve, time_to_maturity, 'linear', 'extrap');


%fprintf('The yield for a time to maturity of %.1f years is %.7f%%.\n', time_to_maturity, yield_at_maturity);
% Define the matrix to store the output variable
outputMatrix = [];

% Print the column headers
fprintf('Time to Maturity\tYield\n');

%For loop that calculates the yield for 
for i = 0:1/12:5
    % Calculate the output variable
    yield_at_maturity = interp1(bond_maturities, NS_Yield_Curve, i, 'linear', 'extrap');

     fprintf('%d/%d\t\t%.6f\n', round(i*12), 12, yield_at_maturity);
end



