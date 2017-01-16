
%% THIS SCRIPT EXPLORES HOW friction ANGLE DISTRIBUTIONS VARY FOR
% FOR DIFFERENT GRAIN SIZE DISTRIBUTIONS. THE CODE WAS MOTIVATED
% BY WIBERG AND SMITH, 1985, KIRCHNER ET AL., 1990 AND BUFFINGTON ET 
% AL., 1992.
    
    % Close and clear items from cache
    close all; clc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP ONE: SPECIFY THE LENGTH OF THE VIRTUAL STREAMBED AND THE 
    % number of grains to randomly drop along the virtual streambed for 
    % friction angle measurement. Do this interactively
    
    prompt = {'Enter the length of the virtual streambed in number of particles (1000 - 5000):',...
        'Enter the number of grains to use for friction angle measurement (5000 - 10000):'};
    dlg_title = 'Friction angle measurement virtual streambed and particle number details';
    num_lines = 1;
    defaultans = {'2500','5000'};
    ILanswer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    % This variable is the length of the virtual streambed as measured in
    % number of particles
    VBed = str2double(ILanswer{1});
    % We need an additional term of +1 the virtual streambed length
    VBedPlus = VBed+1;
    % This variable specifies the number of particles to drop on the
    % virtual streambed for friction angle measurements
    VRandom = str2double(ILanswer{2});    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP TWO: BUILD GRAIN SIZE DISTRIBUTION FOR ANALYSIS
        
        % Call a function to build the grain size distribution. The
        % function builds the distribution based on specifying the mean
        % and length of the distribution measured from the mean outwards.
        [ D50,D84,GSDuse,GS_Dg,GS_Sigma,GSDlength,Grain_mmAvg ]...
            = GSDGeneration_NORMv1( );
                   
            % Now use the results from the specified grain size dist. to
            % figure out some details needed for the calculations below.

            % Specify the minimum grain size in psi units
            PsiMIN = GSDuse(1,2);
            % Specify the maximum grain size in psi units
            PsiMAX = GSDuse(end,2);
            % Find the difference between the max and min psi values
            PsiDIFF = PsiMAX - PsiMIN;
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP THREE: BUILD THE VIRTUAL STREAMBED WITH THE GRAIN SIZE DIST.
            
        %% STEP 3a: POPULATE A LINE SEGMENT WITH VBed GRAINS.  
        % Grain size is determined by a random number, and grain
        % coordinates are specified by their respective diameters. 
        % This is basically a three step process.

            % STEP 3aa: generate a set of random numbers. Shuffle indicates
            % that the random number generator will generate a different 
            % sequence of numbers each time rand is called.

                rng('shuffle');

                % Populate the R vector with VBedPlus random numbers. The 
                % R1 vector requires VBedPlus grains so that a grain 
                % placement address exists just before the last grain b/c 
                % I assume that grains will come to rest at the edge 
                % between two adjacent grains. R1 is used to determine a 
                % random grain size.
                R1 = rand(VBedPlus,1);
                
            %% STEP 3ab: COMPUTE RANDOM GRAIN SIZE FROM RANDOM NO. SEQUENCE
            % Use the random number array to specify random grain sizes to 
            % populate the computational array. Grain size is in mm. 

                % R1GSet is the elevation of the R1GS grain top in mm.
                R1GSet = zeros(VBedPlus,1);

                % Compute R1Psi - R1Psi stands for random grain size in psi
                % units.
                R1Psi = PsiMIN + (PsiDIFF .* R1);
                % Convert the random grain sizes from psi units to mm. 
                % R1GS stands for random grain size in mm.
                R1GS = 2 .^ R1Psi;

            %% STEP 3ac: POPULATE THE LINE SEGMENT WITH RANDOM GRAIN SIZES.
            % Pre-populate the location and elevation vecotrs. centx is the 
            % spatial coordinate of the grain center along the line  
            % segment. Center coordinate is in mm.

                centx = zeros(VBed,1);
                % edgex is the spatial coordiante of far grain edge along  
                % the line segment. Edge coordinate is in mm.
                edgex = zeros(VBed,1);
                % e0 is the elevation of each grain along the grain size 
                % vector, set arbitrarily to 1 so that the center of each  
                % grain lies in the same horizontal plane. Elevation is 
                % in mm.
                e0 = zeros(VBedPlus,1);

                MovingStaCumm = zeros(VBedPlus,VBedPlus);
                MovingStaCummRev = zeros(VBedPlus,VBedPlus);
                constant = zeros(VBedPlus,1);
                
        %% STEP 3b: RUN THROUGH AND BUILD PARAMETER VALUES FOR THE LINE 
        % segment grain locations, from 1 to end of the grain size vector.
               
            for j = 1:VBedPlus

               if j == 1
                   edgex(j) = 0 + R1GS(j);
                   centx(j) = 0 + (R1GS(j) / 2);
                   e0(j) = 1;
                   R1GSet(j) = e0(j) + (R1GS(j) / 2);

               else
                   edgex(j) = edgex(j-1) + R1GS(j);
                   centx(j) = edgex(j-1) + (R1GS(j) / 2);
                   e0(j) = 1;
                   R1GSet(j) = e0(j) + (R1GS(j) / 2);

               end

            end
            clear j
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP FOUR: BUILD FORWARD AND BACKWARD MATRICES OF CUM. STATIONING    
    % This is needed to compute average elevation around any grain based
    % Kirchner et al., 1990 because I will need to search around a grain
    % for a distance of D84 to compute the average grain elevation.
        
        %% STEP 4a: BUILD THE MOVING WINDOW MATRIX OF CUM. POSITION
        % I will read the edgex vector to a new term to keep naming
        % convention consistent through the next several steps. 
        % StaCumm is the cummulative stationing of the random grain size
        % vector.  StaCum is in millimeters.  This term will be used to
        % help compute the average grain top elevation for use in Kirchner
        % et al, 1990.
        StaCumm = edgex;

        % Now compute the moving window matrix of cummulative position
        for j = 1:VBedPlus

           % This if statement creates a matrix of values which provides a
           % moving window of cummulative stationing for the R1GS vector.
           % Negative values in any vector of the matrix will be ignored in
           % the average elevation calculations.
           if j == 1

               % Need this term to create the matrix
               constant(j) = StaCumm(j);
               % This operation subtracts the constant value from each 
               % element of the MovingStaCumm vector.               
               MovingStaCumm(j,:) = StaCumm - constant(j);

           else

               % Need this term to create the matrix
               constant(j) = StaCumm(j);
               % This operation subtracts the constant value from each 
               % element of the MovingStaCumm vector.               
               MovingStaCumm(j,:) = StaCumm - constant(j);

           end  

        end
        clear j
        
        % STEP 4b: BUILD REVERSE MOVING WINDOW MATRIX OF CUM. POSITION
        % - do operation in reverse order.
        
        for j = VBedPlus:-1:1

           if j == VBedPlus
               edgexrev(j) = 0 + R1GS(j);

           else
               edgexrev(j) = edgexrev(j+1) + R1GS(j);

           end

        end
        
        % STEP 4b: BUILD GRAIN ELEVATION MATRICES AND REVERSE MOVING WINDOW 
        % matrix of cum. position - do operation in reverse order.
        
        % I will read the edgexrev vector to a new term to keep naming
        % convenction consistent through the next several steps.
        StaCummRev = edgexrev;
        
        for j = VBedPlus:-1:1

            % This if statement creates a matrix of values which provides a
            % moving window of cummulative stationing for the R1GS vector.
            % Negative values in any vector of the matrix will be ignored
            % in the average elevation calculations.
            if j == VBedPlus
               edgexrev(j) = 0 + R1GS(j);
               % Need this term to create the matrix
               constant(j) = StaCummRev(j);
               % This operation subtracts the constant value from each 
               % element of the MovingStaCumm vector.               
               MovingStaCummRev(j,:) = StaCummRev - constant(j);

           else
               edgexrev(j) = edgexrev(j+1) + R1GS(j);
               % Need this term to create the matrix
               constant(j) = StaCummRev(j);
               % This operation subtracts the constant value from each 
               % element of the MovingStaCumm vector.               
               MovingStaCummRev(j,:) = StaCummRev - constant(j);

           end

        end
        clear j

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP FIVE: COMBINE THE FORWARD AND BACKWARD LOOKING MATRICES INTO 
    % one matrix for later calculations.  I will replace negative values in
    % the forward matrix with positve values from the backward matrix. This
    % is done with logical indexing; I replace the negative components of 
    % the forward looking matrix with positve values of the backward
    % looking matrix.

        MovingStaCumm(MovingStaCumm < 0) = MovingStaCummRev(MovingStaCummRev > 0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP SIX: CREATE A VECTOR OF AVERAGE BED ELEVATION BASED ON A 
    % search neighborhood equivalent to the D84.

        % AvgElevR1GS is the average elevation of the R1 grain size vector
        % based on a search neighborhood of D84. Units are in millimeters.
        AvgElevR1GS = zeros(VBedPlus,1);

        for j = 1:VBedPlus

            [~,c] = find(MovingStaCumm(j,:) < (D84));

            AvgElevR1GS(j) = sum(R1GSet(c)) / length(c);                 

        end
        clear j

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP SEVEN: PLACE GRAINS ALONG LINE SEGMENT & COMPUTE friction ANG.
    % Using the array of VBed random grains placed along the line segment,
    % randomly place grains of random size and compute the friction angle for
    % each occurrence. This is a three step process.

        % STEP 7a: generate two sets of random numbers. Shuffle indicates
        % that the random number generator will generate a different  
        % sequence of numbers each time rand is called.

            rng('shuffle');

            % R2 is used to determine a random grain size.
            R2 = rand(VRandom,1);
            % R3 is used to determine random locations along the grain size 
            % array populated with grains specified by R1. R3 is generated
            % with the random integer generator with replacement, meaning  
            % the same location can be generated more than once.  Note that
            % there are N-1 locations along the array specified above 
            % because it is not possible to compute a friction angle for 
            % the last grain in the sequence because there is no grain down
            % array from it.
            R3 = randi([1,VBed],[VRandom,1]);

        % STEP 7b: locate the R2 random grain sizes along the R1 array
        % according to the random location specified by R3.  Then compute 
        % the R3 grain elevation and finally the friction angle.  Write the
        % friction angle to a parameter and save.

            % R2Psi stands for random grain size in psi units.
            R2Psi = PsiMIN + (PsiDIFF .* R2);
            % Convert the random grain sizes from psi units to mm. 
            % R2GS stands for the randomly placed grain, size expressed in
            % mm.
            R2GS = 2 .^ R2Psi;

        % STEP 7c: map grain sizes from the R1 vector to the R3 location
        % vector for use in looping below and to determine where the placed 
        % grain can physically fit. To determine if the placed grain can
        % physically fit I need to store grain sizes beyond the two 
        % bounding grains of the placement address because their geometry
        % may be incompatible with the size of the placed grain.  In that 
        % case I need to look down array to the next grain size, and 
        % perhaps beyond to determine when the grain may physically fit.

            % R3i term is the R1GS grain size from the R1 vector mapped to
            % the R3 location. This mapping is done by linear indexing with
            % R3. Units are in mm.
            R3GSi = R1GS(R3);
            R3GSiet = (R3GSi ./ 2) + e0(1);
            % R3ii term is the R1GS+1 grain size from the R1 vector mapped
            % to the R3 location. This mapping is done by linear indexing 
            % with R3. Units are in mm.
            R3GSii = R1GS(R3+1);
            % Run through a loop to map the grain sizes upstream of j in
            % order to compute the exposure term properly. the notation of 
            % putting i before the parameter name is to indiacte one value
            % up vector from the i position.
            iR3GSi = zeros(VRandom,1);

            for j = 1:VRandom;

                if R3(j) > 1

                    iR3GSi(j) = R1GS(R3(j)-1);   

                end

            end
            clear j

        % STEP 7d: now roll through and compute friction angles and several 
        % other terms which are needed for evaluating the critical shear 
        % stress function of Wiberg and Smith, 1987; Kirchner et al, 1990.

            dl = zeros(VRandom,1);
            % dl is the length between the
            hyp = zeros(VRandom,1);
            % theta is the computed friction angle in degrees
            theta = zeros(VRandom,1);
            % R2GSe is the elevation of the R2 grain center placed at loc.
            % R3. R2GSet is the elevation of the R2 grain top.
            R2GSec = zeros(VRandom,1);
            R2GSet = zeros(VRandom,1);
            % only use the fist value in the e1 vector because e1 is 
            % constant. if this changes must deal with e1.
            E1(1) = e0(1);
            % exposure is the grain exposure for use in Kirchner1990.
            exposure = zeros(VRandom,1);
            % exposureelevdiff is the relative elevation difference between
            % the random R2 grain top and the upgrain R3i grain top
            exposureelevdiff = zeros(VRandom,1);
            % RelGS is the ratio of the R2 grain size to the D50 of the
            % mixture
            RelGS = zeros(VRandom,1);
            % Define a size ratio condition for use in computing the 
            % friction angles. This is set arbitrarily high so that we can
            % compute frictiona angles for a wide assortment of grain
            % gemoetries. Play with the threshold to examine dependencies.
            SRatio_Threshold = 10000;
            % Compute two grain size ratios used below
            USDS_SRatio = R3GSi ./ R3GSii;
            DSUS_SRatio = R3GSi ./ R3GSii;

            for j = 1:VRandom

                % This if statement tests whether the upstream grain is
                % 1.5x the next downstream grain, and whether the random
                % grain to be placed is larger than the diameters of the 
                % u/s and d/s grains along the R vector.
                if USDS_SRatio(j) >= SRatio_Threshold...
                        || DSUS_SRatio(j) >= SRatio_Threshold

                    if R2GS(j) < ((R3GSi(j) + R3GSii(j)) / 2)

                        % dl computation assumes center of the overlying
                        % grain falls precisely at the edge location  
                        % between the R3i and R3ii grains.
                        dl(j) = (R3GSii(j,1) / 2);

                        % the hypotenuse is computed as the simple sum of 
                        % the radius's for the overlying and R3ii grains.
                        hyp(j) = (R2GS(j) / 2) + (R3GSii(j,1) / 2);

                        % the angle theta is the friction angle.
                        theta(j) = asind(dl(j) / hyp(j));

                    else
                        
                        if USDS_SRatio(j) >= SRatio_Threshold
                        
                            theta(j) = 1;
                            
                        elseif DSUS_SRatio(j) >= SRatio_Threshold
                        
                            theta(j) = 89;
                            
                        end

                    end

                else

                    % dl computation assumes the center of the overlying
                    % grain falls precisely at the edge location between 
                    % the R3i and R3ii grains. units are in mm.
                    dl(j) = (R3GSii(j,1) / 2);

                    % the hypotenuse is computed as the simple sum of the 
                    % radius's for the overlying and the R3ii grains.
                    % units are in mm.
                    hyp(j) = (R2GS(j) / 2) + (R3GSii(j,1) / 2);

                    % the angle theta is the friction angle.
                    theta(j) = asind(dl(j) / hyp(j));

                    % R2GSec is the elevation of the R2 grain center placed 
                    % at location R3 or R3new. R2GSet is the elevation of 
                    % the R2 grain top. units are in mm.
                    R2GSec(j) = ((cosd(theta(j)) * (R2GS(j) / 2))...
                        + (sind(theta(j)) * (R3GSii(j) / 2))) + E1;
                    R2GSet(j) = R2GSec(j) + (R2GS(j) / 2);   

                    % Compute the relative grain size for R2GS
                    RelGS(j) = R2GS(j) ./ (D50);

                    % Compute R2GS grain exposure for use in Kirchner1990v3.
                    % units are in mm.  Exposure equals the top elevation
                    % of the random placed grain minus the top elevation of
                    % bed grain j.  This is a simplification of reality.  
                    % I should additionally search up grain from grain j to
                    % assess if grains extend above random placed grain.  
                    % I will do that in the future.

                    exposureelevdiff(j) = R2GSet(j) - R3GSiet(j);

                    if exposureelevdiff(j) > 0

                        exposure(j) = R2GSet(j) - R3GSiet(j);

                    else

                        exposure(j) = 0;

                    end

                end

            end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP EIGHT: CALL KIRCHNER FUNCTION TO COMPUTE CRITICAL DIMENSIONLESS
    % shear stress distributions from the friction angle distributions

    % Remember that some terms are in units of mm. Convert to meters inside
    % the Kirchner funtion.

        % Set z equal to R2GSet. units are in mm.
        z = R2GSet;
        % Call the modified Kirchner function.
        [taucrit,taucritless,R2GSeval,ycrit,Scrit,Ustarcrit,UStarCrit]...
            = Kirchner1990v3(R2GS,theta,exposure,z,R3GSiet,VRandom,D84);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP NINE: CALL THE MOBILITY FUNCTION TO COMPUTE PDFS OF MOBILITY
    % based on 0.5 PSI class intervals.
        
        % Call the mobility PDF function. Output of this function describes
        % the percentile values of the dimensionless Shields Stress for
        % percentile classes: 10 20 30 40 50 60 70 80 and 90.
        [TCMPercentileM,Percentile]...
            = MobilityPDFMOB_v1(GSDlength,taucrit,taucritless,PsiMIN,PsiMAX,R2GSeval,VRandom);

        % Build a matrix to store the critical dimensionless stress
        % percentile results by grain size.  The results are stored
        % as row = grain size and column = critical dimensionless
        % stress percential class
        TCMPercentileMatrix = TCMPercentileM;
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP TEN: SAVE SOME VARIALBES TO FILE FOR USE ELSEWHERE
    
    DepthCrit = ycrit;
    SlopeCrit = Scrit;
    UshearCrit = Ustarcrit;
    UShearCrit = UStarCrit;
    ShearCrit = taucrit;
    Shields = taucritless;
    ShieldsPercMatrix = TCMPercentileMatrix;
            
    save('Shields.mat','Shields');
    save('Stress.mat','ShearCrit');
    save('DepthCrit.mat','DepthCrit');
    save('SlopeCrit.mat','SlopeCrit');
    save('UshearCrit.mat','Ustarcrit');
    save('ShieldsPercentiles','ShieldsPercMatrix');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FINAL STEP: CREATE SOME FIGURES
    
    % Tile the average grain size vector to match the matrix of
    % dimensionless Shields stresses from the MobilityPDF function
    [Per,Pec] = size(Percentile);
    Grain_mmAvg_Mat = repmat(Grain_mmAvg',1,Pec);
    [Pr,Pc] = size(Grain_mmAvg_Mat);
    Percentile_Mat = repmat(Percentile,Pr,1);
      
    
    % Define the colormap
    cmap = jet(length(taucritless));
    colormap(cmap);
    
    % Define a plotting vector to make log scale colorbars
    cc = [0.0010 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1.0];
    
    %% FIRST FIGURE
    h2 = figure(2);
    for j = 1:Pc
        scatter(Grain_mmAvg_Mat(:,j),ShieldsPercMatrix(:,j),[],Percentile_Mat(:,j),'filled','MarkerEdgeColor','k')
        hold on
    end
    box on
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    dim = [0.76 0.7 0.3 0.3];
    str = {'Percentile Class'};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor',[1 1 1]);
    xlabel('Average Grain Diameter (mm)','fontsize',12,'fontweight','b');
    ylabel('\tau_c_r*','fontsize',12,'fontweight','b');
    title('Simulated Distribution of Shields Stress');
    colorbar('FontSize',8,'YTickLabel',Percentile);
    
    FAModelFig2 = ('ShieldsPercentile.png');
    % Export a file and specify dimensions, etc
    set(gcf, 'Units','centimeters')
    set(gcf,'PaperType','B5','PaperPositionMode','auto') 
    print(h2,'-dpng',FAModelFig2,'-r600');
    
    %% NEXT FIGURE   
    h3 = figure(3);
    data1 = log10(RelGS);
    data2 = log10(theta);
    data3 = log10(taucritless);
    scatter(data1,data2,[],data3,'filled','MarkerEdgeColor','k')
    title('Simulated Distribution of Friction Angle');
    box on
    dim = [0.73 0.7 0.3 0.3];
    str = {'Critical Shields Stress'};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor',[1 1 1]);
    xlabel('log_{10} D_i / D_5_0','fontsize',12,'fontweight','b');
    ylabel('log_{10} Friction angle \theta','fontsize',12,'fontweight','b');
    caxis(log10([cc(1) cc(end)]));
    colorbar('FontSize',8,'YTick',log10(cc),'YTickLabel',cc);
    
    FAModelFig3 = ('ShieldsDistribution.png');
    % Export a file and specify dimensions, etc.
    set(gcf, 'Units','centimeters')
    set(gcf,'PaperType','B5','PaperPositionMode','auto') 
    print(h3,'-dpng',FAModelFig3,'-r600');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Last step
    msgbox('Probabilistic mobility script completed successfully');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  



