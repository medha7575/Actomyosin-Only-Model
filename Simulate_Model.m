%% SIMULATION STARTS HERE
% This file is used to run the simulation after the initialization is
% completed. It first calculates force on individual
% nodes of actomyosin material and then moves moves nodes based on the force applied, thus implementing Equation 1-4.
% Finally it generates output in form of images to visualize and
% study material behavior in the given setup.

NodeDensity = zeros(4000,2);
for mcs=1:3500
    mcs
    
    
    NewActin=0;
    
    F_am=zeros(Mc_NodeCount,2);  % net Actomyosin force
    F_r_Nb= zeros(Mc_NodeCount,MaxNeighAllowed); %  Actomyosin force in x direction by ith neighbor
    F_c_Nb= zeros(Mc_NodeCount,MaxNeighAllowed);  %  actomyosin force in y direction by ith neighbor
    F_r_NbSp= zeros(Mc_NodeCount,MaxNeighAllowed); % Spring force in x direction by ith neighbor
    F_c_NbSp= zeros(Mc_NodeCount,MaxNeighAllowed);  % Spring force in y direction by ith neighbor
    F_Spring=zeros(Mc_NodeCount,2);     % Net PRing force by all nbeighbors
    F_myoT=zeros(Mc_NodeCount,2);
    
    
    %% ----------Myosin forces calulation------
    
    RandNodes=randperm(Mc_NodeCount);
    for i = 1: Mc_NodeCount
        NodeID= RandNodes(i);
        %   %   Calculating forces
        if Mc_Node(NodeID,MOVEORNOT)==1
            %             disp('found a fixed node, so net force is zero')
            F_MyoT(NodeID,ROW)= 0;
            F_MyoT(NodeID,COL)= 0;
            % disp('NodeID=')NodeID ;
        else
            % Force calculation of myosin  starts here....
         
            for i = 1: Mc_Node(NodeID,NO_OF_NEIGHBOR)  % for all neighbors
                
                NbID = neighborhood(NodeID,i);  % get Neighbor iD
                
                tn_c_Nb= Mc_Node(NbID,COL)- Mc_Node(NodeID,COL);  % OUTWARD VECTOR since pulling from neighbor will be directed out
                
                tn_r_Nb= Mc_Node(NbID,ROW)- Mc_Node(NodeID,ROW);
                
                tn_Nb = 0.001+sqrt((tn_r_Nb)^2 + (tn_c_Nb)^2);   % distance between node and neighbor (adding .001 to avaoid numerical instaibility.)
                
                
                
                %% If there is tension in the spring
                %   ----Implementing Equation-2 ---------
                F_NbSp=0;
                if tn_Nb>tn_rest
                    
                    % % Spring force calculation
                    
                    F_NbSp = KSp*(tn_Nb- tn_rest); % force generated due to spring = contstant times spring tension
                    
                    F_r_NbSp(NodeID,i)= F_NbSp*(tn_r_Nb/tn_Nb); % outward force felt by the node due to tension
                    F_c_NbSp(NodeID,i)= F_NbSp*(tn_c_Nb/tn_Nb);% This is equal to magnitude*unit vector
                    
                end
                
                
                F_Spring(NodeID,ROW)=  F_Spring(NodeID,ROW) + F_r_NbSp(NodeID,i);%adding forces by all the neighbors
                F_Spring(NodeID,COL)=  F_Spring(NodeID,COL)+ F_c_NbSp(NodeID,i);
                
                % % Inter-myosin force calculation
                   %   ----Implementing Equation-1 ---------
            % % calculate pulling forces from  myosin neighbors
                F_Nb=0;
                AM_Neighbor= Mc_Node(NbID,MYO_CONC);  %fetch the myosin activity  on node and neighbor 
                AM_BaseNode= Mc_Node(NodeID,MYO_CONC);
                
                F_Nb =(K*AM_Neighbor*AM_BaseNode); % force magnitude
                
                F_r_Nb(NodeID,i)=F_Nb*(tn_r_Nb/tn_Nb);  % Actomyosin force on the node due to ith nieghbor
                F_c_Nb(NodeID,i)=F_Nb*(tn_c_Nb/tn_Nb); % magnitude * unit vector 
                
                
                F_am(NodeID,ROW)= F_am(NodeID,ROW)+ F_r_Nb(NodeID,i); %%adding forces by all the neighbors
                F_am(NodeID,COL)= F_am(NodeID,COL)+ F_c_Nb(NodeID,i); %
                
                
                
            end
            F_myoT(NodeID,ROW)= F_am(NodeID,ROW)+F_Spring(NodeID,ROW);
            F_myoT(NodeID,COL)= F_am(NodeID,COL)+F_Spring(NodeID,COL);
            
        end
        
    end
    
    %% ----Following lines of code apply the forces calculated above to move the nodes
    %----- Implementing equation-4------
    
    Combo = zeros(bound,bound); %% Resetting the lattice
    
    
    RandNodes=randperm(Mc_NodeCount);
    for i = 1: Mc_NodeCount
        NodeID= RandNodes(i);
   
        colNo= Mc_Node(NodeID,COL);
        rowNo= Mc_Node(NodeID,ROW);
        
        
        % % rounding to 2 didgit to control spatial resolution
        NewcolNo = round(colNo+ F_myoT(NodeID,COL),2);
        NewrowNo = round(rowNo+ F_myoT(NodeID,ROW),2);
        
        %%Next each node checks if the space is available or already
        %%saturated by maximum no of nodes 
        noOfNodesAtNewPosition = NodeCountPerPixel((round(NewrowNo)),(round(NewcolNo)));
        
        if noOfNodesAtNewPosition == Max_myosin_nodes_per_pixels
            NewcolNo= colNo;
            NewrowNo= rowNo;
        end
        
        
%%----- Update the variables with new position to implement the move---
        Mc_NodeNo(round(Mc_Node(NodeID,ROW)),round(Mc_Node(NodeID,COL)))= -1;   %% reset the NodeNo in order to create it again afresh
        
        Mc_Node(NodeID,ROW)= NewrowNo; %
        Mc_Node(NodeID,COL)= NewcolNo;
        
        Mc_NodeNo(round(NewrowNo),round(NewcolNo))=NodeID;  
        
        Combo(round(NewrowNo),round(NewcolNo))= Mc_Node(NodeID,COMBOVALUE); % refresh combo
        
        NodeCountPerPixel((round(NewrowNo)),(round(NewcolNo))) = NodeCountPerPixel((round(NewrowNo)),(round(NewcolNo))) + 1;
        NodeCountPerPixel((round(rowNo)), (round(colNo))) = NodeCountPerPixel((round(rowNo)),(round(colNo))) - 1 ;
        
        
    end
    %%%%%%%%%%%%%%%%%%------END--------%%%%%%%%%%%
    
    %% ------Visualize the output ------
    
    if (mcs)==2
        close all
        cd (absoluteFolderPath)
        Actinfoldername=strcat('Combo');
        mkdir(Actinfoldername);
        cd (Actinfoldername)
        filename=strcat('Actin_',num2str(mcs));
        fig = figure;
        C=uint8(Combo);
        imshow(C)
        axis square
        print(fig,filename,'-dpng');
        
        cd ..
        cd ..
        
        %
    end
    
    if mod(mcs,100)==0
        close all
        cd (absoluteFolderPath)
        Actinfoldername=strcat('Combo');
        mkdir(Actinfoldername);
        cd (Actinfoldername)
        filename=strcat('Actin_',num2str(mcs));
        fig = figure;
        %                 imagesc(Combo)
        %                 colormap('hot')
        C=uint8(Combo);
        imshow(C)
        %                      imshow(C, 'Colormap', jet(255))
        axis square
        print(fig,filename,'-dpng');
        
        cd ..
        cd ..
        
        %
    end
    
    %
end