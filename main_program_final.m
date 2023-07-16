clear 
clc

rng default
load('LMM_model_calibrated');

students_decide_to_do_optional_additional_section = false;
num_simulations = 1;

for(k=1:num_simulations)


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Section:(Total 5 marks) You have to create a strategy to decide how
% much money is invested in the bonds (Assets) portfolio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%AUM = 8bn for APAC
t                               = 1;
AUM(t)                          = 8000000000; % AUM = Assets under management
fixed_amount_extract            = 0.01*AUM(t); % monthly fixed ammount that investors extract every month

%We choose the amount of money invested in bonds further on in the sim


size_historical_term_structure = length(RateSpecG2); %From 2005 to 2013 


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Section:(Total 5 marks) You need to decide which periods are ran by 
% your simulation, and you have to "create" a strategy for
% including default risk (credit risk spread)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We add in credit spreads in section 4

N = 3;      % num of years

%%%========================= TO BE DONE! ======================
% Random start selection -- You need to use a strategy for the selection of the
% first day between 10-May-2005 and 06-Jul-2011  (2 Marks)

start = RateSpecG2{1}.Settle;

% Load historical yield curve data from CSV file
opts = detectImportOptions('yield-curve-rates-1990-2021.csv');
opts1 = detectImportOptions('2023-treasury-rates-7apr.csv');
opts.VariableNamingRule = 'preserve';
opts1.VariableNamingRule = 'preserve';
opts = setvaropts(opts, 'Date', 'InputFormat', 'MM/dd/yy');
opts1 = setvaropts(opts1, 'Date', 'InputFormat', 'MM/dd/yyyy');

historical_yield_curves = readtable('yield-curve-rates-1990-2021.csv', opts); 

% Filter out dates outside of the range of 2005 to 2011
historical_yield_curves = historical_yield_curves(year(historical_yield_curves.Date) >= 2005 & year(historical_yield_curves.Date) <= 2011, :);

todays_yield_curve = readtable('2023-treasury-rates-7apr.csv', opts1); 

% Extract the interest rate values for 2 years and 10 years:
two_years = todays_yield_curve{:, 6}; % 6th column represents '2 Yr' maturity
ten_years = todays_yield_curve{:, 10}; % 11th column represents '10 Yr' maturity
today_slope = ten_years - two_years; % calculate the slope

% Calculate the slope of each historical yield curve
historical_slopes = historical_yield_curves{'2005-01-03' <= historical_yield_curves.Date & historical_yield_curves.Date <= '2011-12-30', '10 Yr'} - historical_yield_curves{'2005-01-03' <= historical_yield_curves.Date & historical_yield_curves.Date <= '2011-12-30', '2 Yr'}; 

% Calculate the difference between each historical slope and today's slope 
differences = abs(historical_slopes - today_slope); 

% Find the historical yield curve with the smallest difference from today's yield curve 
[min_difference, min_index] = min(differences); 

% Extract the date with the closest yield curve
closest_dates = historical_yield_curves.Date(find(differences == min_difference));



% differences minimised on 15th and 16th Nov 2006
% 15th is best as rates nearer todays. 4.80% 2yr and 4.61% 10yr
myDate = datenum(closest_dates(2)-1);

current_day_ir_simulation = myDate - start;
%%%============================================================


list_days_simulation = busdays(RateSpecG2{1}.Settle+current_day_ir_simulation+1,RateSpecG2{1}.Settle+current_day_ir_simulation+365*N); % from next day to the next 3 years 


initial_day             = list_days_simulation(1);
current_month           = month(initial_day);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Section:(Total 5 marks) You have to create a strategy to decide how
% much money is invested in the bonds (Assets) portfolio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Moved first section here to load in variables that I use to select
%starting cash amount
months_remaining_simulation  = 36; 

%Creating vector last_bus_date which contains all the dates we pay out
%fixed amount
day_one = initial_day;
for i = 1:months_remaining_simulation
        last_bus_date(i) = busdate(eomdate(day_one), 1);
        % Add one month to the current date
        next_month = addtodate(day_one, 1, 'month');
        day_one = next_month;
end

%Calculating vector dif which is days between current date and future
%payment dates
dif = last_bus_date - initial_day;

%rate is a vector of the annualised treasry rates for the time to each payment
rate =  getZeroRates(RateSpecG2{current_day_ir_simulation},RateSpecG2{current_day_ir_simulation}.Settle + dif);
months_left_vector = 1:months_remaining_simulation;

%Calculating the present value of each payment
value = fixed_amount_extract./(1+rate'/12).^months_left_vector;
PV_Liability                    = sum(value);

%Only investing in bonds the amount that we need to cover PV_liability
portion_to_maintain_in_cash     = 1 - (PV_Liability/AUM(t));
%%%============================================================
actual_cash_amount              = AUM(t)*portion_to_maintain_in_cash;
actual_investment_amount        = AUM(t)*(1-portion_to_maintain_in_cash); % cash + bonds investment = AUM 
amount_to_subtract_interest_up  = 100000;
amount_to_add_interest_down     = 100000;

%%%========================= TO BE DONE! ======================
% The 10 year interest rate (r(0,10)) was selected for that random day (previous sub-section). Is that
% enough? HINT: What can be done to simulate Credit Spreads? (Remember the
% variable ``RateSpecG2'' contains the interest rate term structure of
% treasuries securities, that has different risk than the bonds that you
% might be using. (3 Marks)
R10Y_Note = getZeroRates(RateSpecG2{current_day_ir_simulation},RateSpecG2{current_day_ir_simulation}.Settle + 10*365);
previous_rate_10Y_note   = R10Y_Note;  

%getting the day diference to calculate exactly one month previosu
%this method is more effective since it takes into account leap years and
%months with different numbeer days (28,29,30,31)
%finding the index in RateSpecG2 for the one month previous day and getting
%the 10Y note from then
one_month_previous_day = datenum(datestr(datetime(datestr(initial_day)) - calmonths(1)));
days_difference = initial_day - one_month_previous_day;
one_month_previous_day_ir_simulation = current_day_ir_simulation - days_difference;
previous_rate_10Y_note = getZeroRates(RateSpecG2{one_month_previous_day_ir_simulation}, RateSpecG2{one_month_previous_day_ir_simulation}.Settle + 10*365);

%%%============================================================



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third Section:(Total 10 marks) You have to adjust the suggested ``fictitious'' bonds with the ``real'' selected bonds, and 
% you have to decide the amount to be invested in every bond the first day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%========================= TO BE DONE! ======================
% Adjust the suggested ``fictitious'' bonds with the ``real'' selected bonds (TO BE DONE!) (10 Marks)
%Read in real bonds from excel file
data = readtable('bonds.xlsx');


% Extract the different variables from the excel file as vectors
P = table2array(data(1, 2:end));
coupon = table2array(data(2, 2:end));
face_value = table2array(data(3, 2:end));
bond_maturity = table2array(data(4, 2:end));
coupons_per_year = table2array(data(5, 2:end));
bond_basis = table2array(data(6, 2:end));


bond_maturity_dates = zeros(1,length(bond_maturity)); % Preallocate the vector for the maturity dates
for i = 1:length(bond_maturity)
    bond_maturity_dates(i) = list_days_simulation(1) + bond_maturity(i)*365;
end



bond_units = [];


%immunising the portfolio at the start of the simulation
Settle = myDate;

%calculating different variables
for i = 1:length(coupon)
        CouponRate  = coupon(i);
        Maturity    = bond_maturity_dates(i);
        ZeroDates   = [(initial_day):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(i);
        Period      = coupons_per_year(i);
        Basis       = bond_basis(i);
        start_yield(i)        =  bndyield(P(i),CouponRate, Settle, Maturity,'Period',Period,'Basis',Basis, 'Face', FaceValue); %getting the yield of the different bonds at the start of simulation
 end

for(i=1:length(coupon))
    %improved duration calculation by adding in additional
    %variables and using a more accurate r value
    [ModDuration(i),Duration(i),PerDuration(i)] = bnddury(start_yield(i),coupon(i),Settle,bond_maturity_dates(i),coupons_per_year(i),bond_basis(i),'Face', face_value(i));
end
%adding in convexity to calculation
for(i=1:length(coupon))
    [Convexity(i),PerConvexity(i)] = bndconvy(start_yield(i),coupon(i),Settle,bond_maturity_dates(i),coupons_per_year(i),bond_basis(i),'Face', face_value(i));
end

num_total_bonds         = length(coupon);
num_combination_bonds   = 10 ; % Here we select the number of bonds of our portfolio
index_bond = nchoosek([1:num_total_bonds],num_combination_bonds); % This creates the list of ALL the possible combinations of bonds to immunise the portfolio

Maturity_Portfolio  = months_remaining_simulation/12;
for(j = 1:size(index_bond,1)) % WE test ALL the combinations. You have to select one from this for loop
        
    % Matrix for solutions

   
    A = [ones(1,num_combination_bonds); ...
        1/PV_Liability*[Duration(index_bond(j,:))]; ...
        1/PV_Liability*[Convexity(index_bond(j,:))] ];
    
    b = [PV_Liability;...
        Maturity_Portfolio;...
        Maturity_Portfolio^2];
    
    W(j,:) = A\b;
    
    Units(j,:) = round(W(j,:)./(P(index_bond(j,:))));
    Value_portfolio(j,:) = P(index_bond(j,:)).*Units(j,:);
    TE(t,j) =  sum(Value_portfolio(j,:)) - PV_Liability;
end

%selecting strategy with minimum tracking error
[~, immunisation_strategy_selected] = min(abs(TE(t,:)));

if(~any(isnan(Units(immunisation_strategy_selected,:))))
    list_of_bonds_to_use            = index_bond(immunisation_strategy_selected,:);
    list_of_bonds_not_to_use        = setdiff(1:num_total_bonds,index_bond(immunisation_strategy_selected,:));
    bond_units(t,list_of_bonds_to_use) =  Units(immunisation_strategy_selected,:);
    bond_units(t,list_of_bonds_not_to_use) =  0;

end


amount_to_invest        = sum(bond_units(t,:).*P);
actual_cash_amount      = actual_cash_amount + (actual_investment_amount-amount_to_invest); % We add the remaining cash from the rounding of the bond units to the cash account
cash(t)                 = actual_cash_amount;
actual_investment_amount= amount_to_invest; 
portfolio_bonds_current_value(t) = actual_investment_amount;
%%%============================================================

disp(['Day ' datestr(list_days_simulation(t)) ' -- AUM: $'  num2str(AUM(t)) ' -- Cash amount: ' num2str(cash(t)) ' -- Bonds Portfolio Value:' num2str(portfolio_bonds_current_value(t))]);



%%
t = 2; % We start the simulation from the second day


for(current_date = list_days_simulation(2:end)')
%for(current_date = list_days_simulation(1:30)')
    if(month(current_date)~= current_month)
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fourth Section: (Total 10 marks) Cash inflows/outflows (They happen at the
        % beginning of the month
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % We are starting a new month. We need to settle cash outflows (pension fund retirments) and
        % cash inflows (pension fund deposits), and update the value of the portfolio
        current_month = month(current_date);
        months_remaining_simulation = months_remaining_simulation-1;% we subtract one month from the remaining months of the simulation
        last_bus_date = last_bus_date(2:end); %remove the most recent payment date from the vector as it will have been paid out below
        
        %%%========================= TO BE DONE! ======================
        % The 10 year interest rate (r(0,10)) was selected for that random day (previous sub-section). Is that
        % enough? (HINT: What can be done to simulate Credit Spreads? Remember the
        % variable ``RateSpecG2'' contains the interest rate term structure of
        % treasuries securities, that has different risk than the bonds that you
        % might be using. (4 Marks)
        treasury = getZeroRates(RateSpecG2{current_day_ir_simulation},RateSpecG2{current_day_ir_simulation}.Settle + bond_maturity_dates - current_date);

        spread = (treasury' - r); %calculating the credit spreads between our bonds and the treasury bonds
        
        investment_sum        = sum(abs(bond_units(t-1,:)).*P');
        
        weighting = investment_sum/sum(investment_sum); %finding weighting of portfolio
        
        portfolio_spread = weighting * spread'; %creating a weighted average for our credit spreads


        R10Y_Note = getZeroRates(RateSpecG2{current_day_ir_simulation},RateSpecG2{current_day_ir_simulation}.Settle + 10*365);
        current_rate_10Y_note = R10Y_Note + portfolio_spread; %adding credit spread to ten year treasury note



        Settle      = current_date;
        for i = 1:length(coupon)
                CouponRate  = coupon(i);
                Maturity    = bond_maturity_dates(i);
                ZeroDates   = [(current_date+1):365:Maturity];
                ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
                FaceValue   = face_value(i);
                Period      = coupons_per_year(i);
                Basis       = bond_basis(i);
                [Price_t(i)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
                r(i)        =  bndyield(Price_t(i),CouponRate, Settle, Maturity,'Period',Period,'Basis',Basis, 'Face', FaceValue);
         end
         value_of_bonds = dot(bond_units(t-1,:),Price_t);
        %%%============================================================
        
        if(current_rate_10Y_note-previous_rate_10Y_note>0.001) % Increase of of more than 10bp
            
            if(actual_cash_amount-amount_to_subtract_interest_up>0) % We need to check we can extract such amount of cash
                
                
                
                actual_cash_amount = actual_cash_amount-(current_rate_10Y_note-previous_rate_10Y_note)*amount_to_subtract_interest_up;
            else
                %%%========================= TO BE DONE! ======================
                % You have to sell bonds to get cash for paying the
                % retirements of the fund (TO BE DONE!) (3 Marks)
                %%%============================================================
                cash_required = amount_to_subtract_interest_up - actual_cash_amount;

               
                %here we are finding the ratio of bonds needed to be
                %sold. we then sell a ratio of those bonds, show by the
                %subtraction in bond units and increase in actual cash
                %amount
                ratio_of_bonds_needed_to_sell = cash_required / value_of_bonds;
                bond_units = bond_units - ratio_of_bonds_needed_to_sell * bond_units;
                actual_cash_amount = actual_cash_amount + ratio_of_bonds_needed_to_sell * value_of_bonds;
                
                % Then, you can subtract the cash
                actual_cash_amount = actual_cash_amount-amount_to_subtract_interest_up;
            end
                
        elseif(current_rate_10Y_note-previous_rate_10Y_note<-0.001) % Decrease of more than 10bp
           
            actual_cash_amount = actual_cash_amount+(previous_rate_10Y_note-current_rate_10Y_note)*amount_to_add_interest_down; % We don't need to check cash amount when it is added
            
        else
            % Nothing changes in terms of the fixed cash inflows/outflows
            
        end
        if(actual_cash_amount-fixed_amount_extract>0)% We need to check we can extract such amount of cash
            actual_cash_amount = actual_cash_amount - fixed_amount_extract;
        else
            %%%========================= TO BE DONE! ======================
            % You have to sell bonds to get cash for paying the
            % retirements of the fund (TO BE DONE!) (3 Marks)
            %%%============================================================

            cash_required = fixed_amount_extract - actual_cash_amount;

            
           
            %here we are finding the ratio of bonds needed to be
            %sold. we then sell a ratio of those bonds, show by the
            %subtraction in bond units and increase in actual cash
            %amount
            ratio_of_bonds_needed_to_sell = cash_required / value_of_bonds;
            bond_units = bond_units - ratio_of_bonds_needed_to_sell * bond_units;
            actual_cash_amount = actual_cash_amount + ratio_of_bonds_needed_to_sell * value_of_bonds;
            
            % Then, you can subtract the cash
            actual_cash_amount = actual_cash_amount-fixed_amount_extract; 
        end
        previous_rate_10Y_note  = current_rate_10Y_note;% update  current rate for the next step. Only changes once per month
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fifth Section (Total 60 marks): We change the composition of the portfolio: You
        % re-immunise the portfolio at the beginning of every month
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%========================= TO BE DONE! ======================        
        % Adjust the suggested ``fictitious'' bonds with the ``real'' selected bonds (TO BE DONE!) (5 Marks)
       
       
        %%%============================================================
        
        %%%========================= TO BE DONE! ======================        
        % Improve the calculation of the duration (10 Marks)
        Settle = current_date;
        for(i=1:length(coupon))
            %improved duration calculation by adding in additional
            %variables and using a more accurate r value
            [ModDuration(i),Duration(i),PerDuration(i)] = bnddury(r(i),coupon(i),Settle,bond_maturity_dates(i),coupons_per_year(i),bond_basis(i),[],[],[],[],[],face_value(i));
        end
        %%%============================================================
        
        %%%========================= TO BE DONE! ======================        
        % Include convexity in the immunisation 
        for(i=1:length(coupon))
            [Convexity(i),PerConvexity(i)] = bndconvy(r(i),coupon(i),Settle,bond_maturity_dates(i),coupons_per_year(i),bond_basis(i),[],[],[],[],[],face_value(i));
        end

        num_total_bonds         = length(coupon);
        num_combination_bonds   = 10 ; % Here we select the number of bonds of our portfolio
        index_bond = nchoosek([1:num_total_bonds],num_combination_bonds); % This creates the list of ALL the possible combinations of bonds to immunise the portfolio
        
        %%%========================= TO BE DONE! ======================
        % Improve the estimation of the PV of Liability (8 Marks)        
        % Imporved Pv of liability by adjusting the r value for each
        % iteration of the fixed payments we are obliged to
        dif = last_bus_date - current_date;
        if length(dif) ~= 0
            rate =  getZeroRates(RateSpecG2{current_day_ir_simulation},RateSpecG2{current_day_ir_simulation}.Settle + dif);
            months_left_vector = 1:months_remaining_simulation;
            value = fixed_amount_extract./(1+rate'/12).^months_left_vector;
            PV_Liability                    = sum(value);
        else 
            %PV_liability goes to 0 during last month when we have made all
            %payments
            PV_Liability = 0;
        end


        Maturity_Portfolio  = months_remaining_simulation/12;
        
        for(j = 1:size(index_bond,1)) % WE test ALL the combinations. You have to select one from this for loop
            
            % Matrix for solutions

           
            A = [ones(1,num_combination_bonds); ...
                1/PV_Liability*[Duration(index_bond(j,:))]; ...
                1/PV_Liability*[Convexity(index_bond(j,:))] ];
            
            b = [PV_Liability;...
                Maturity_Portfolio;...
                Maturity_Portfolio^2];
            
            W(j,:) = A\b;
            
            Units(j,:) = round(W(j,:)./(Price_t(index_bond(j,:))));
            Value_portfolio(j,:) = Price_t(index_bond(j,:)).*Units(j,:);
            TE(t,j) =  sum(Value_portfolio(j,:)) - PV_Liability;
        end
        %%%============================================================

        %%%========================= TO BE DONE! ======================        
        % In here you have to select the immunised portfolio.
        % We suggested the second, but this is not optimal.
        %%% EXTREMELY IMPORTANT: YOU HAVE TO EXPLAIN WHY YOU SELECTED SUCH
        %%% STRATEGY!!!!! (Include such explanation in the report)  (10 Marks)
        %%% ==========================================================
        %Select the strategy with the minimum tracking error
        [~, immunisation_strategy_selected] = min(abs(TE(t,:)));
        
        if(~any(isnan(Units(immunisation_strategy_selected,:))))
            list_of_bonds_to_use            = index_bond(immunisation_strategy_selected,:);
            list_of_bonds_not_to_use        = setdiff(1:num_total_bonds,index_bond(immunisation_strategy_selected,:));
            bond_units(t,list_of_bonds_to_use) =  Units(immunisation_strategy_selected,:);
            bond_units(t,list_of_bonds_not_to_use) =  0;
            
            % Change in number of bonds
            cash_change = sum((bond_units(t,:)-bond_units(t-1,:)).*Price_t);
            actual_cash_amount = actual_cash_amount - cash_change;
        else
            % No change is made. It's not possible
            bond_units(t,:) = bond_units(t-1,:);
        end
        %%%============================================================

        if(students_decide_to_do_optional_additional_section)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Seventh Section: (Total 10 Marks): The current bond prices
            % provide a precise scenario of the interest rate term
            % structure. Nevertheless, the interest rate term structure is
            % dynamic. You might use the code below to improve the
            % immunisation code by including Monte Carlo simulations of the
            % possible future scenarios.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%========================= OPTIONAL! ======================
            
            nPeriods = 10;
            nTrials  = 10000;
            [LMMZeroRatesSimPaths, LMMForwardRatesSimPaths] = LMM{current_day_ir_simulation}.simTermStructs(nPeriods+1,'nTrials',nTrials,'antithetic',true);
            %%%============================================================
        end
    else
        %% ? (TO BE DONE?)
        
        
        bond_units(t,:) = bond_units(t-1,:);
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sixth Section: (Total 10 Marks): Updates in the value of the portfolio 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%========================= TO BE DONE! ======================
    % Adjust with the selected bonds (TO BE DONE!) (2 Marks)
    Settle      = current_date;

    for i = 1:length(coupon)
        CouponRate  = coupon(i);
        Maturity    = bond_maturity_dates(i);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(i);
        Period      = coupons_per_year(i);
        Basis       = bond_basis(i);
        [Price_t(i)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(i) = bndyield(Price_t(i),CouponRate, Settle, Maturity,'Period',Period,'Basis',Basis, 'Face', FaceValue);
    end

   
    
    
    %%%============================================================
    
    
    % Current Valuation:
    portfolio_bonds_current_value(t)= sum(bond_units(t,:).*Price_t); % Where P_t is the price of bonds today (TO BE DONE!)
    cash(t)                         = actual_cash_amount;
    
    %%%========================= TO BE DONE! ======================
    % Improve the estimation of the PV of Liability (8 Marks)
    Liability                       = months_remaining_simulation*fixed_amount_extract; % total fixed extracts per month remainings;
    Maturity_Portfolio              = months_remaining_simulation/12;
    m                               = 12; % compounds per year  

    %Pv of liability is improved the same as earlier
    dif = last_bus_date - current_date;

   if length(dif) ~= 0
        rate =  getZeroRates(RateSpecG2{current_day_ir_simulation},RateSpecG2{current_day_ir_simulation}.Settle + dif);
        months_left_vector = 1:months_remaining_simulation;
        value = fixed_amount_extract./(1+rate'/12).^months_left_vector;
        PV_Liability                    = sum(value);
    else 
        PV_Liability = 0;
    end
    %PV_Liability                    =  Liability/(1+r_average/m).^(Maturity_Portfolio*m);    
    %%%============================================================
    actual_TE(t,1)                         =  portfolio_bonds_current_value(t) - PV_Liability;
    
    
    AUM(t) = portfolio_bonds_current_value(t) + cash(t);
    
    disp(['Day ' datestr(list_days_simulation(t)) ' -- AUM: $'  num2str(AUM(t)) ' -- Cash amount: ' num2str(cash(t)) ' -- Bonds Portfolio Value:' num2str(portfolio_bonds_current_value(t))]);
    %disp(bond_units)
    % We increase the simulation day counters  by 1
    current_day_ir_simulation = current_day_ir_simulation+1;
    t = t+1;
end


TE_simulation{k} = actual_TE;

end


