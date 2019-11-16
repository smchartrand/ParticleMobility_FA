function [ D50,D84,GSDuse,GS_Dg,GS_Sigma,GSDlength,Grain_mmAvg ]...
    = GSDGeneration_NORMv1( )

    %% THIS CODE CREATES A NORMALLY DISTRIBUTED GRAIN SIZE DISTRIBUTION
    % SPECIFYING A GEOMTRIC MEAN GRAIN SIZE AND THE STANDARD DEVIATION.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 1: CREATE NORMAL DISTRIBUTION
    
    prompt = {'Enter the geometric mean grain size in millimeters (2.0 - 32.0):',...
        'Enter the standard deviation of the distribution (0.25 - 2.0):'};
    dlg_title = 'Friction angle measurement virtual streambed and particle number details';
    num_lines = 1;
    defaultans = {'7.3','0.68'};
    ILanswer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    % Pass the distribution mean to a variable
    GS_Dg = str2double(ILanswer{1});
    % Convert the geometric mean value (in millimeters) to Psi units
    GS_DgPsi = log2(GS_Dg);
    %GS_DgPsi = log2(GS_Dg) / log(2);
    % Pass the standard deviation input to a variable
    GS_Sigma = str2double(ILanswer{2});
    % Specify the limits x of the normal distribution - this is equivalent
    % to defining the Psi-based limits of the grain size distribution pdf
    Upper_GSLimit = GS_DgPsi + (4 * GS_Sigma);
    Lower_GSLimit = GS_DgPsi - (4 * GS_Sigma);
    ND_Limits = (Lower_GSLimit:0.01:Upper_GSLimit);
    % Convert the limits to grain size in millimeters
    GS_inMM = 2.^ND_Limits;
    % Make the distribution based on the information just provided
    Distro = cdf('Normal',ND_Limits,GS_DgPsi,GS_Sigma);
    % Plot the distribution for inspection
    % Specify a plot style
    plotStyle = {'b-'};
    % Make the plot
    h1 = figure(1);
    plot(GS_inMM,Distro,plotStyle{1});
    % Make x-axis log scale
    set(gca,'xscale','log')
    ax = gca;
    % x-axis limits
    ax.XLim = [10^(-1) 10^(2.5)];
    % Turn grid on
    grid on;
    % Specify axis labels
    ylabel('Percent finer'),xlabel('Grainsize (mm)');
    title('Simulated Grain Size Distribution');
    
    FAModelFig1 = ('Model_GSDistribution2.png');
    % Export a file and specify dimensions, etc.
    set(gcf, 'Units','centimeters')
    set(gcf,'PaperType','B5','PaperPositionMode','auto') 
    print(h1,'-dpng',FAModelFig1,'-r600');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 2: COMPUTE THE GRAIN SIZE STATISTICS NEEDED FOR THE FRICTION 
    % ANGLE MEASUREMENTS. THIS REQUIRES THREE COLUMN VECTORS WITH GRAIN SIZE 
    % (MM), GRAIN SIZE (PSI) AND THE ASSOCIATED CUMULATIVE DISTRIBUTION VALUE
    
    % Specify Psi Scale
    PsiLowLimit = -4.0;
    PsiUppLimit = 12.0;
    PsiScale = PsiLowLimit:0.5:PsiUppLimit;
    
    % Specify the PsiRange vector for analysis
    PsiLow = PsiScale(find(PsiScale<Lower_GSLimit,1,'last'));
    PsiUpp = PsiScale(find(PsiScale>Upper_GSLimit,1,'first'));
    %PsiRange = (PsiLowLimit:0.5:PsiUppLimit);
    PsiRange = (PsiLow:0.5:PsiUpp);
    [~,PsiN] = size(PsiRange);
    if  PsiN == 1
        error('Error. Grain size distribution has only one grain size')
    end
    % Find indices which line up with half-scale psi values over the range
    % specified above, and write the associated grain sizes to a new
    % vector, and then write the associated cumulative distribution
    % probability density values to a new vector
    % Index variable
    Index(1,length(PsiRange)) = zeros;
    % Grain size distribution in mm
    GS_Distro_MM(1,length(PsiRange)) = zeros;
    % Grain size cumulative probability distribution
    GS_Distro_CD(1,length(PsiRange)) = zeros;
    % Enter loop
    for j = 1:length(PsiRange)
        
        if j < length(PsiRange)
            
           Index(j) = find(ND_Limits > PsiRange(j),1,'first');
            
        else
            
           Index(j) = find(ND_Limits < PsiRange(j),1,'last'); 
            
        end
        GS_Distro_MM(j) = GS_inMM(1,Index(j));
        GS_Distro_CD(j) = Distro(1,Index(j));
    
    end
    
    % Identify the D50 grain size    
    D50 = GS_inMM(mean(find(Distro < 0.5,1,'last'),find(Distro > 0.5,1,'first')));
    % Identify the D84 grain size 
    D84 = GS_inMM(mean(find(Distro < 0.84,1,'last'),find(Distro > 0.84,1,'first')));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP THREE: PACKAGE THE GRAIN SIZE AND CUMULATIVE DISTRIBUTION 
    % VECOTRS INTO A MATRIX 
    
    % Create the matrix and flip the data so it lines up with the friction
    % angle calculations
    GSD = [GS_Distro_MM; PsiRange; GS_Distro_CD]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP FOUR:IDENTIFY THE GRAIN SIZE (MM AND PSI) AND THE CUMULATIVE 
    % DISTRIBUION TO USE FOR THE FRICTION ANGLE CALCULATIONS
            
    % First identify the address of the 'last' ~ 0 value in the cumm
    % distribution of the computed grain size distribution and then
    % the 'first' 1 in the distribution.  This eliminates size
    % classes which do not account for mass within the distribution
    g0 = find(GSD(:,3)>0.000001,1,'first');
    g100 = find(GSD(:,3)>0.9999,1,'first');
    % Changed this to operate on the array rather than row or
    % column for each spatial node.
    GSDuse = GSD(g0:g100,:);
    
    % Compute a few other statistics before closing this out
    GSDlength = length(GSDuse) - 1;
    % Set the vector index parameter - the number of grain classes    
    j = 1:GSDlength;
    % Compute characteristic grain size of each size class in the bed.
    % First compute the average Psi class for the grain size distro
    Grain_PsiAvg(j) = ((GSDuse(j,2) + GSDuse(j+1,2)) * 0.5);
    % Now compute the size in mm.
    Grain_mmAvg(j) = 2 .^ Grain_PsiAvg(j);
    clear j
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
