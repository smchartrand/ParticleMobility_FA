function [ taucrit,taucritless,R2GSeval,ycrit,Scrit,Ustarcrit,UStarCrit ] = Kirchner1990v3( R2GS,theta,exposure,z,R3GSiet,VRandom,D84 )

%% THIS FUNCTION COMPUTES CRITICAL DIMENSIONLESS SHEAR STRESS BASED ON THE 
% work of Kirchner et al., 1990, which stems from Wiberg and Smith, 1987.
% Critical dimensionless shear stress is computed principally from a 
% measure of a grain friction angle.  V2 of the code reflects a
% modification of the original Eq. 12 ecause it was found to not work for
% cases of poorly sorted sediment.  This resulted in re-formulation of the
% integral in Eq. 12.

    %% DEFINE A FEW CONSTANTS RATHER THAN PASSING THEM TO THE FUNTION
    % Density of water - kg per cubic meter
    rhow = 1000;
    
    % Density of sediment - kg per cubic meter
    rhos = 2650;

    % Gravitational acceleration - meters per square second
    g = 9.80665;
    
    % Manning-Strickler alpha parameter - dimensionaless
    a = 8.1;
    
    % Pi
    Pi = pi();
    
    % Drag coefficient - dimensionless (source Wiberg and Smith, 1985)
    Cd = 0.40;
    
    % Lift coefficient - dimensionless (source Wiberg and Smith, 1985)
    Cl = 0.20;
    
    % Von Karmen constant - dimensionless
    K = 0.407;
    
    % Kinematic viscosity of water - meters / sec.
    KinVis = 1.14 * 10^-6;
    
    % Pre-populate the relevant terms.
    fz = zeros(VRandom,1);
    eD = zeros(VRandom,1);
    feD = zeros(VRandom,1);
    Term1 = zeros(VRandom,1);
    Term2 = zeros(VRandom,1);
    Term3 = zeros(VRandom,1);
    Term4 = zeros(VRandom,1);
    Term5 = zeros(VRandom,1);
    taucrit = zeros(VRandom,1);
    taucritless = zeros(VRandom,1);
    ReGrain = zeros(VRandom,1);
    evaluate = zeros(VRandom,1);
    R2GSeval = zeros(VRandom,1);
    ycrit = zeros(VRandom,1);
    Scrit = zeros(VRandom,1);
    Ustarcrit = zeros(VRandom,1);
    
    % Convert projection and exposure to meters
    exposurem = exposure ./ 1000;
    R2GSm = R2GS ./ 1000;
    R3GSietm = R3GSiet ./ 1000;
    zm = z ./ 1000;
    
    Ksx = 2 * D84;
    % For 42 lps, average width of 0.58 meters.
    qwpass = 0.0707;
    % Average initial bed slope
    S = 0.018;
    
    % Compute the normal depth approximation according to Parker e-book
    ynorm = ((Ksx^0.33*qwpass^2)/(a^2*g*S))^0.30;
        
    %% COMPUTE THE CRITICAL SHEAR FUNCTION IN PIECES TO AVOID MISTAKES.
    % I have broken the function into five terms.  The third term 
    % involves the integral from p-e to p, which I solve using the 
    % trapezoidal rule.  The other terms are algebraic.
    
    % Role through one for loop to compute the trapeziodal x stations, the 
    % f(z) term, and the trapezoidal function terms for use in the 
    % approximation below.
    for i = 1:VRandom
        
        % Compute the f() functions. Careful with units.

        fz(i) = log((zm(i) + R3GSietm(i)) ./ R3GSietm(i));

        eD(i) = ((((exposurem(i) - R2GSm(i)) + R3GSietm(i)) ./ R3GSietm(i)));
        
        if exposurem(i) <= 0.005 || eD(i) < 0
            
            feD(i) = 0;
            
        else

            feD(i) = log(eD(i));
            
        end
        
        % Compute term 1 of Kirchner et al., 1990 Eq. 12
        Term1(i) = (rhos - rhow) .* g .* ((Pi .* (R2GSm(i).^3)) ./ 6);

        % Compute term 2 of Kirchner et al., 1990 Eq. 12
        Term2(i) = 0.5 .* (tand(theta(i)).^-1) .* Cd .* ((K.^2)^-1); 

        % Compute term 3 of Kirchner et al., 1990 Revised form of Eq. 12        
        Term3(i) = Pi .* (R2GSm(i) ./ 2) .* exposurem(i) .* fz(i).^2;

        % Compute term 4 of Kirchner et al., 1990 Eq. 12
        Term4(i) = Pi .* Cl .* ((K.^2)^-1) .* (R2GSm(i) ./ 2) .* exposurem(i); 

        % Compute term 5 of Kirchner et al., 1990 Eq. 12
        Term5(i) = (fz(i).^2);% - (feD(i).^2); 

        % Compute the critical dimensional and dimensionless stress.  Also
        % compute the Reynolds grain number, the Stokes number and the
        % Settling velocity
        evaluate(i) = Term1(i) .* ((Term2(i) .* Term3(i)) + (Term4(i) .* Term5(i))).^-1;
        
        if exposurem(i) == 0 || evaluate(i) <= 0
            
            taucrit(i) = NaN;
            taucritless(i) = NaN;
            ReGrain(i) = NaN;
            R2GSeval(i) = NaN;
            ycrit(i) = NaN;
            Scrit(i) = NaN;
            Ustarcrit(i) = NaN;
            
        else                
                
            taucrit(i) = Term1(i) .* ((Term2(i) .* Term3(i)) + (Term4(i) .* Term5(i))).^-1;
            taucritless(i) = taucrit(i) ./ (( rhos - rhow) .* g .* R2GSm(i));
            ReGrain(i) = ((sqrt(taucrit(i)./rhow)) .* R2GSm(i)) ./ KinVis;
            % Compute the critical depth
            ycrit(i) = (taucritless(i) .* ( rhos - rhow) .* g .* R2GSm(i)) ./ (rhow .* g .* 0.018);
            % Compute the critical slope
            Scrit(i) = (taucritless(i) .* ( rhos - rhow) .* g .* R2GSm(i)) ./ (rhow .* g .* ynorm);
            Ustarcrit(i) = (g * ycrit(i) * Scrit(i))^0.5;
            R2GSeval(i) = R2GSm(i) .* 1000;
            
        end
        
    end
    
    if i == VRandom
       
       % This is the average shear velocity for all grain size and grain positions 
       MeanSCrit = nanmedian(Scrit);
       MeanYCrit = nanmedian(ycrit);
       MeanTau = nanmedian(taucritless);
       UStarCrit = (g * MeanYCrit * MeanSCrit) ^ 0.5;
        
    end

end

