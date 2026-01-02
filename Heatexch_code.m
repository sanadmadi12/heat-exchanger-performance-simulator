%% Simulator/ main scripts
disp(" ");
disp("====================================");
disp("   HEAT EXCHANGER INPUT PANEL");
disp("====================================");
disp(" ");
disp("Enter 1 for a single-tube heat exchanger.")
disp("Enter 2 for a multipass shell-and-tube heat exchanger.")
disp("Enter 3 for a crossflow heat exchanger (single pass, both fluids unmixed)")
r=input("Enter your choice (1–3): ");
if r==1
    [Q,T1h,T2h,T1c,T2c,regime1,regime2,y,flow,Tln]=simulator();
    p = NaN;                % not used in r==1
    printResults(Q, T1h, T2h, T1c, T2c,r, y, p, flow,regime1, regime2,Tln, NaN);
    
%% Case 2
elseif r==2
    disp("Choose the shell-and-tube configuration:")
    disp("1) One shell pass (number of tube passes must be a multiple of 2)")
    disp("2) Multiple shell passes (shell passes must be a multiple of 2; tube passes a multiple of 4)")   
    p=input("Enter your choice (1–2): ");
    if p==1 %1 shell
        [Q,T1h,T2h,T1c,T2c,regime1,regime2,y,flow,Tln]=simulator();
        [P,R]=RandP(T1h,T2h,T1c,T2c);
        sqR=sqrt(R^2 + 1);
        F=(sqR * log((1 - P)/ (1 - R*P)))/((R - 1)*log((2-P*(R+1 -sqR))/(2-P*(R+1+sqR))));
        
      
        if ~isreal(F)||isnan(F)||isinf(F)||F<=0||F>1
            F=0.75;
        end
        Qf = Q*F;
        if Qf<0
            Qf=Qf*-1;
        end
        printResults(Qf, T1h, T2h, T1c, T2c,r, y, p, flow,regime1, regime2,Tln, F);
    
    elseif p==2 %2 shell

        [Q,T1h,T2h,T1c,T2c,regime1,regime2,y,flow,Tln]=simulator();
        [P,R]=RandP(T1h,T2h,T1c,T2c);

        X = (1 - P*R) / (1 - P);

        if abs(R - 1) < 1e-6
        % Special case if R is exactly 1
            P_prime = P / (2 - P);
        else
        % Standard formula for P_prime
            P_prime = (1 - X^(1/2)) / (R - X^(1/2));
        end
        sqR = sqrt(R^2 + 1);

        if abs(R - 1) < 1e-6
            % Special case for R = 1
            sq2 = sqrt(2);
            F = (sq2 * P_prime / (1 - P_prime)) / log((2 - P_prime*(2 - sq2)) / (2 - P_prime*(2 + sq2)));
        else
            % Standard formula using P_prime
            F = (sqR * log((1 - P_prime) / (1 - R*P_prime))) /((R - 1) * log((2 - P_prime*(R + 1 - sqR)) / (2 - P_prime*(R + 1 + sqR))));
        end

        if ~isreal(F)||isnan(F)||isinf(F)||F<=0||F>1
            F=0.75;
        end
        Qf = Q*F;
        if Qf<0
            Qf=Qf*-1;
        end

        printResults(Qf, T1h, T2h, T1c, T2c,r, y, p, flow,regime1, regime2,Tln, F);
        
        
    else 
        error("Invalid input. Please enter a valid option.")
    end
    

elseif r==3
    [Q,T1h,T2h,T1c,T2c,regime1,regime2,y,flow,Tln]=simulator();
    [P,R]=RandP(T1h,T2h,T1c,T2c);
    pfun=@(n) (1 - exp((n^0.22 / R) * (exp(-R * n^0.78) - 1)))-P;
    ntu=fzero(pfun,[0.001,50]);
    if abs(R - 1) < 1e-4
        % Use the limit formula for R = 1
        F = P / (ntu * (1 - P));
    else
        % Use the standard formula
        F = log((1 - R*P) / (1 - P)) / (ntu * (R - 1));
    end
    if ~isreal(F) || isnan(F) || isinf(F) || F<=0 || F>1
        F = 0.75;
    end
    Qf=Q*F;

    if Qf<0
        Qf=Qf*-1;
    end

    p = NaN; % not used in r==3
    printResults(Qf,T1h,T2h,T1c,T2c,r,y,p,flow, regime1,regime2, Tln,F);

else 
    error("Invalid input")
end
%% functions
%1
function[R]=Rvalue(hi,ho,rfi,rfo,Ai,Ao,Di,Dto,kwall,L)
    R1=1/(hi*Ai);
    R2=1/(ho*Ao);
    R3=rfi/Ai;
    R4=rfo/Ao;
    R5 = log(Dto/Di)/(2*pi*kwall*L);
    R=R1+R2+R3+R4+R5;
end 

%2
function[Re]=Revalue(V,x,rho,D)
    v=x/rho;
    Re=V*D/v;
end

%3
function [Nu]= Nuvalue(Re,Pr,s)
    if Re<2300
        if s==1
            disp("Tube-side flow is laminar. Choose a boundary condition: ")
            disp("1) Constant heat flux")
            disp("2) Constant surface temperature")
            type=input("Enter your choice (1–2): ");
            if type==1
                Nu=4.36;
            elseif type==2
                Nu=3.66;
            else
                error("Invalid input. Please enter a valid option.")
            end
        else
            disp("Shell-side flow is laminar. Choose a boundary condition: ")
            disp("1) Constant heat flux")
            disp("2) Constant surface temperature")
            type=input("Enter your choice (1–2): ");
            if type==1
                Nu=4.36;
            elseif type==2
                Nu=3.66;
            else
                error("Invalid input. Please enter a valid option.")
            end
        end
    elseif Re >= 2300 && Re <= 4000
        phi=(Re-2300)/(4000-2300);
        Nulam=3.66;
        f = 0.3164*Re^(-0.25);
        Nuturb = ((f/8) * (Re - 1000) * Pr) /(1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1));
        Nu=(1-phi)*Nulam+phi*Nuturb;
    else
        f = 0.3164*Re^(-0.25);
        Nu = ((f/8) * (Re - 1000) * Pr) /(1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1));
    end
end
% 
function regime = flowRegime(Re)
    
    if Re < 2300
        regime = "Laminar";
    elseif Re <= 4000
        regime = "Transitional";
    else
        regime = "Turbulent";
    end
end

%4
function[h]=hvalue(Nu,k,D)
    h=(k/D)*Nu;
end

%5
function[Q]=QknownT(Ui,Ai,T1,T2,T3,T4)
    DT1=T1-T2;
    DT2=T3-T4;
    if abs(DT1-DT2) < 1e-9
        Tln=abs(DT1);   % LMTD limit as DT1 -> DT2
    else
        Tln=abs(abs(DT1-DT2)/log(abs(DT1/DT2)));
    end
    
    Q=Ui*Ai*Tln;
end

%6
function[m,c]=MandC(m1,m2)
    disp("Select the mass flow rate for the hot fluid:")
    disp(["1) ṁ =" num2str(m1) " kg/s"])
    disp(["2) ṁ =" num2str(m2) " kg/s"])
    b=input("Enter your choice (1–2): ");
    if b==1
        m=m1;
    elseif b==2
        m=m2;
    else
        error("Invalid input. Please enter 1 or 2.")
    end
    c=input("Enter the specific heat capacity, c (J/kg·K): ")
end

%7
function[T1h,T2h,T1c,T2c,Q]=missingTCounterQ(m,c,U,A)
    disp("Which temperature is unknown?")
    disp("1) T_h,out (hot outlet)")
    disp("2) T_c,out (cold outlet)")
    disp("3) T_h,in (hot inlet)")
    disp("4) T_c,in (cold inlet)")
    miss = input("Enter your choice (1–4): ");
    if miss==1
        T1h=input("Enter hot inlet temperature, T_h,in: ")
        T2h=NaN;
        T1c=input("Enter cold inlet temperature, T_c,in: ")
        T2c=input("Enter cold outlet temperature, T_c,out: ")
        if T1c>T2c
            error("Invalid temperatures: cold outlet must be greater than cold inlet.")
        end        
        F=@(T2h) (m*c*(T1h-T2h))- U*A*(((T1h - T2c) - (T2h - T1c)) / log((T1h - T2c) / (T2h - T1c)));

        a = T1c+1e-6;   
        b = T1h-1e-6;

        T2h=fzero(F,[a b]);
        disp("Solved hot outlet temperature, T_h,out = " + num2str(T2h,"%.6g"));

        Q=m*c*(T1h-T2h);
        disp("Heat transfer rate, Q = " + num2str(Q,"%.6g") + " W")
      

    elseif miss==2
        T1h=input("Enter hot inlet temperature, T_h,in: ")
        T2h=input("Enter hot outlet temperature, T_h,out: ")
        T1c=input("Enter cold inlet temperature, T_c,in: ")
        T2c=NaN;
        if T1h<T2h
            error("Invalid temperatures: hot outlet cannot be greater than hot inlet.")
        end

        Q=m*c*(T1h-T2h);
        
        F=@(T2c) Q- U*A*(((T1h - T2c) - (T2h - T1c)) / log((T1h - T2c) / (T2h - T1c)));

        a = T1c + 1e-6;
        b = T1h - 1e-6;

        T2c=fzero(F,[a b]);
        


    elseif miss==3
        T1h=NaN;
        T2h=input("Enter hot outlet temperature, T_h,out: ")
        T1c=input("Enter cold inlet temperature, T_c,in: ")
        T2c=input("Enter cold outlet temperature, T_c,out: ")
        if T1c>T2c
            error("Invalid temperatures: cold outlet must be greater than cold inlet.")
        end
        F=@(T1h) (m*c*(T1h-T2h))- U*A*(((T1h - T2c) - (T2h - T1c)) / log((T1h - T2c) / (T2h - T1c)));

        a = max(T2h, Tcc) + 1e-6;
        b = a+200;
        T1h = fzero(F,[a b]);
        

        Q=m*c*(T1h-T2h);
        
    elseif miss==4
        T1h=input("Enter hot inlet temperature, T_h,in: ")
        T2h=input("Enter hot outlet temperature, T_h,out: ")
        T1c=NaN;
        T2c=input("Enter cold outlet temperature, T_c,out: ")
        if T1h<T2h
            error("Invalid temperatures: hot outlet cannot be greater than hot inlet.")
        end
        Q=m*c*(T1h-T2h);
        
        F=@(T1c) Q- U*A*(((T1h - T2c) - (T2h - T1c)) / log((T1h - T2c) / (T2h - T1c)));

        b = min([T1h, T2h, T2c]) - 1e-6;
        a = b - 200;
        T1c = fzero(F,[a b]);
        
    else 
        error("Invalid input. Please enter a numbFer from 1 to 4.")
    end
end

%8
function[T1h,T2h,T1c,T2c,Q]=missingTParallelQ(m,c,U,A)
    disp("Which temperature is unknown?")
    disp("1) T_h,out (hot outlet)")
    disp("2) T_c,out (cold outlet)")
    disp("3) T_h,in (hot inlet)")
    disp("4) T_c,in (cold inlet)")
    miss = input("Enter your choice (1–4): ");
    if miss==1
        T1h=input("Enter hot inlet temperature, T_h,in: ")
        T2h=NaN;
        T1c=input("Enter cold inlet temperature, T_c,in: ")
        T2c=input("Enter cold outlet temperature, T_c,out: ")
        if T1c>T2c
            error("Invalid temperatures: cold outlet must be greater than cold inlet.")
        end
        
        F=@(T2h) (m*c*(T1h-T2h))- U*A*(((T1h - T1c) - (T2h - T2c)) / log((T1h - T1c) / (T2h - T2c)));

        a = T2c + 1e-6;  
        b = T1h - 1e-6;

        T2h = fzero(F,[a b]);
        

        Q=m*c*(T1h-T2h);
        
      

    elseif miss==2
        T1h=input("Enter hot inlet temperature, T_h,in: ")
        T2h=input("Enter hot outlet temperature, T_h,out: ")
        T1c=input("Enter cold inlet temperature, T_c,in: ")
        T2c=NaN;
        if T1h<T2h
            error("Invalid temperatures: hot outlet cannot be greater than hot inlet.")
        end


        Q=m*c*(T1h-T2h);
        
        F=@(T2c) Q- U*A*(((T1h - T1c) - (T2h - T2c)) / log((T1h - T1c) / (T2h - T2c)));

        a = T1c + 1e-6;
        b = T2h - 1e-6;

        T2c = fzero(F,[a b]);
        


    elseif miss==3
        T1h=NaN;
        T2h=input("Enter hot outlet temperature, T_h,out: ")
        T1c=input("Enter cold inlet temperature, T_c,in: ")
        T2c=input("Enter cold outlet temperature, T_c,out: ")
        if T1c>T2c
            error("Invalid temperatures: cold outlet must be greater than cold inlet.")
        end

        F=@(T1h) (m*c*(T1h-T2h))- U*A*(((T1h - T1c) - (T2h - T2c)) / log((T1h - T1c) / (T2h - T2c)));

        a = max(T2h, T1c) + 1e-6;
        b = a+200;
        T1h = fzero(F,[a b]);
        

        Q=m*c*(T1h-T2h);
       

    elseif miss==4
        T1h=input("Enter hot inlet temperature, T_h,in: ")
        T2h=input("Enter hot outlet temperature, T_h,out: ")
        T1c=NaN;
        T2c=input("Enter cold outlet temperature, T_c,out: ")
        if T1h<T2h
            error("Invalid temperatures: hot outlet cannot be greater than hot inlet.")
        end

        Q=m*c*(T1h-T2h);
        
        F=@(T1c) Q- U*A*(((T1h - T1c) - (T2h - T2c)) / log((T1h - T1c) / (T2h - T2c)));

        b = min(T1h, T2c,T2h)-1e-6;
        a =b-200;

        T1c = fzero(F,[a b]);
        
    else 
        error("Invalid input")
    end
end

%9
function[P,R]=RandP(T1,T2,T3,T4)
    R=(T1-T2)/(T4-T3);
    P=(T4-T3)/(T1-T3);
end

%10
function printResults(Q,T1h,T2h,T1c,T2c, r,y,p,flow, regime1,regime2, Tln,F)

    % ---- Safe defaults so nothing crashes ----
    if nargin < 1 || isempty(Q),      Q = NaN; end
    if nargin < 2 || isempty(T1h),    T1h = NaN; end
    if nargin < 3 || isempty(T2h),    T2h = NaN; end
    if nargin < 4 || isempty(T1c),    T1c = NaN; end
    if nargin < 5 || isempty(T2c),    T2c = NaN; end

    if nargin < 6 || isempty(r),      r = NaN; end
    if nargin < 7 || isempty(y),      y = NaN; end %#ok<NASGU>
    if nargin < 8 || isempty(p),      p = NaN; end
    if nargin < 9 || isempty(flow),   flow = NaN; end

    if nargin < 10 || isempty(regime1), regime1 = ""; end
    if nargin < 11 || isempty(regime2), regime2 = ""; end
    if nargin < 12 || isempty(Tln),     Tln = NaN; end
    if nargin < 13 || isempty(F),       F = NaN; end

    regime1 = string(regime1);
    regime2 = string(regime2);

    disp(" ");
    disp("====================================");
    disp("      HEAT EXCHANGER RESULTS");
    disp("====================================");

    % ---------- CONFIGURATION ----------
    flowconfig = "";

    if ~isnan(r)
        if r == 1
            if ~isnan(flow)
                if flow == 1
                    flowconfig = "Single tube – Counterflow";
                elseif flow == 2
                    flowconfig = "Single tube – Parallel flow";
                else
                    flowconfig = "Single tube";
                end
            else
                flowconfig = "Single tube";
            end

        elseif r == 2
            if ~isnan(p)
                if p == 1
                    flowconfig = "Shell-and-tube – 1 shell pass";
                elseif p == 2
                    flowconfig = "Shell-and-tube – 2 shell passes";
                else
                    flowconfig = "Shell-and-tube";
                end
            else
                flowconfig = "Shell-and-tube";
            end

        elseif r == 3
            flowconfig = "Crossflow – both fluids unmixed";
        end
    end

    if strlength(flowconfig) > 0
        disp("Configuration:");
        disp("  " + flowconfig);
    end

    % ---------- HEAT TRANSFER ----------
    if ~isnan(Q)
        disp(" ");
        disp("Heat transfer:");
        disp("  Q = " + num2str(Q,"%.6g") + " W");
    end

    % ---------- TEMPERATURES ----------
    temps = [T1h T2h T1c T2c];
    if all(~isnan(temps))
        disp(" ");
        disp("Temperatures:");
        disp("  Hot side:");
        disp("    Th,in  = " + num2str(T1h,"%.6g"));
        disp("    Th,out = " + num2str(T2h,"%.6g"));
        disp("  Cold side:");
        disp("    Tc,in  = " + num2str(T1c,"%.6g"));
        disp("    Tc,out = " + num2str(T2c,"%.6g"));
    end

    % ---------- FLOW REGIMES ----------
    if strlength(regime1) > 0 || strlength(regime2) > 0
        disp(" ");
        disp("Flow regimes:");
        if strlength(regime1) > 0
            disp("  Side 1: " + regime1);
        end
        if strlength(regime2) > 0
            disp("  Side 2: " + regime2);
        end
    end

    % ---------- LMTD ----------
    if ~isnan(Tln)
        disp(" ");
        disp("Thermal performance:");
        disp("  LMTD = " + num2str(Tln,"%.6g"));
    end

    % ---------- F (ONLY IF PROVIDED) ----------
    if ~isnan(F)
        disp("  F = " + num2str(F,"%.6g"));
    end

    disp("====================================");
    disp(" ");
end


%% sim
function [Q,T1h,T2h,T1c,T2c,regime1,regime2,y,flow,Tln]=simulator()
    regime1 = "";
    regime2 = "";
    Tln = NaN;
    flow = NaN;
    y = NaN;
    %% Heat Exhg properites 
    Di=input("Enter the inner tube diameter, D_i (m): ");
    Do=input("Enter the outer tube diameter, D_o (m): ");
    if Di>=Do
        error("Invalid input: outer diameter must be greater than inner diameter.")
    end
    L=input("Enter the tube length, L (m): ");
    disp("1) Compute Q using known convection coefficients (h)")
    disp("2) Compute Q using fluid properties (solve for h)")
    y=input("Enter your choice (1–2): ");
    if y==1
        hi=input("Enter the inner-side convection coefficient, h_i (W/m²·K): ")
        ho=input("Enter the outer-side convection coefficient, h_o (W/m²·K): ")
        rfi=input("Enter the inner-side fouling factor, R_fi (m²·K/W). Enter 0 if none: ")
        rfo=input("Enter the outer-side fouling factor, R_fo (m²·K/W). Enter 0 if none: ")
        mdot1=input("Enter the tube-side mass flow rate, ṁ_1 (kg/s): ")
        %c1=input("input the scpecific heat capacity for fluid in tube: ")
        mdot2=input("Enter the shell-side mass flow rate, ṁ_2 (kg/s): ")
        %c2=input("input the scpecific heat capacity for fluid in shell: ")
        
        twall=input("Enter the tube wall thickness, t (m). Enter 0 for a negligible wall: ")
        if twall==0
            kwall=1;
        else 
            kwall=input("Enter the wall thermal conductivity, k_wall (W/m·K): ")
        end 
        

        Dto =Di+2*twall;
        Ai=pi*Di*L;
        Ao=pi*Dto*L;

        R=Rvalue(hi,ho,rfi,rfo,Ai,Ao,Di,Dto,kwall,L);
        

        Ui=1/(R*Ai);
        Uo=1/(R*Ao);
        regime1="";
        regime2="";

    elseif y==2
    Aicross=pi*(Di^2)/4;
    Ai=pi*Di*L;

    disp("Enter the following properties for the tube-side fluid (Fluid 1):")
    rho1=input("Density, ρ (kg/m³): ")
    x1=input("Dynamic viscosity, μ (Pa·s): ")
    mdot1=input("Mass flow rate, ṁ (kg/s): ")
    %c1=input("input the scpecific heat capacity: ")
    k1=input("Thermal conductivity, k (W/m·K): ")
    Pr1=input("Prandtl number, Pr: ")
    
    V1=mdot1/(rho1*Aicross);
    
    disp("Enter the following properties for the shell-side fluid (Fluid 2):")
    rho2=input("Density, ρ (kg/m³): ")
    x2=input("Dynamic viscosity, μ (Pa·s): ")
    mdot2=input("Mass flow rate, ṁ (kg/s): ")
    %c2=input("input the scpecific heat capacity: ")
    k2=input("Thermal conductivity, k (W/m·K): ")
    Pr2=input("Prandtl number, Pr: ")

    rfi=input("Enter the inner-side fouling factor, R_fi (m²·K/W). Enter 0 if none: ")
    rfo=input("Enter the outer-side fouling factor, R_fo (m²·K/W). Enter 0 if none: ")

    twall=input("Enter the tube wall thickness, t (m). Enter 0 for a negligible wall: ")
        if twall==0
            kwall=1;
        else 
            kwall=input("Enter the wall thermal conductivity, k_wall (W/m·K): ")
        end
    Dto=Di+2*twall;
    Dh=Do-Dto;

    A2=pi*(Do^2-Dto^2)/4;

    V2=mdot2/(rho2*A2);
    

    Re1=Revalue(V1,x1,rho1,Di);
    regime1=flowRegime(Re1);
    Re2=Revalue(V2,x2,rho2,Dh);
    regime2=flowRegime(Re2);

    Nu1=Nuvalue(Re1,Pr1,1); %tube
    Nu2=Nuvalue(Re2,Pr2,2); %shell

    h1=hvalue(Nu1,k1,Di); %tube
    h2=hvalue(Nu2,k2,Dh); %shell


     
    Ui = 1/((1/h1)+ rfi+ (Di/(2*kwall))*log(Dto/Di)+(Di/Dto)*((1/h2) + rfo )); %tube
    Uo = 1 /((Dto/Di) *( (1/h1) + rfi )+(Dto/(2*kwall)) *log(Dto/Di) +((1/h2) + rfo)); %shell
    end
    %% Calculate Q

    disp("Enter 1 if all four temperatures are known.")
    disp("Enter 2 if three temperatures are known and the fourth is unknown.")
    x=input("Enter your choice (1–2): ");
    if x==1 %%All temps are known
        T1h=input("Enter hot inlet temperature, T_h,in: ")
        T2h=input("Enter hot outlet temperature, T_h,out: ")
        T1c=input("Enter cold inlet temperature, T_c,in: ")
        T2c=input("Enter cold outlet temperature, T_c,out: ")
        if T1h<T2h || T1c>T2c
            error("Invalid temperatures: hot outlet cannot exceed hot inlet, and cold outlet must exceed cold inlet.")
        end
        disp("Enter 1 for counterflow (always select for multipass and crossflow).")
        disp("Enter 2 for parallel flow.")
        flow=input("Enter your choice (1–2): ");
        if flow==1
            Q=QknownT(Ui,Ai,T1h,T2c,T2h,T1c);
        elseif flow==2
            Q=QknownT(Ui,Ai,T1h,T1c,T2h,T2c);
        else 
            error("Invalid input")
        end
    elseif x==2 %3 temps are known

        disp("Enter 1 for counterflow (always select for multipass and crossflow).")
        disp("Enter 2 for parallel flow.")
        flow=input("Enter your choice (1–2): ");
        if flow==1 %counter flow
            [m,c] = MandC(mdot1, mdot2);
            [T1h,T2h,T1c,T2c,Q] = missingTCounterQ(m,c,Ui,Ai);   
        elseif flow==2 %parallel flow
            [m,c] = MandC(mdot1, mdot2);
            [T1h,T2h,T1c,T2c,Q]=missingTParallelQ(m,c,Ui,Ai);
        else 
            error("Invalid input")
        end
        
    else 
        error("Invalid input")
    end
    if Q<0
        Q=Q*-1;
    end 
    
    if flow==1      % counterflow
        DT1=T1h-T2c;
        DT2=T2h-T1c;
    elseif flow==2  % parallel
        DT1=T1h-T1c;
        DT2=T2h-T2c;
    end

    if abs(DT1-DT2)<1e-9
        Tln=abs(DT1);
    else
        Tln=abs(abs(DT1-DT2)/log(abs(DT1/DT2)));
    end    
end


