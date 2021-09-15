%% VARIANT OF SHAMIR's scan (inspired by the SPACE FILLING CURVE)
% BASED ON A CONTROLLED RANDOM WALK
% Author: Kannan Karthik
% Oct 24, 2018


clear all;
close all;
% GRID FORMATION
N = 21; % AN ODD NUMBER
% INDICATES AN N X N grid of points which will be traversed by a GRAPH
C = [(N+1)/2, (N+1)/2]; % CENTRE OF THE GRID FROM WHICH THE RANDOM WALK will ORIGINATE
GRID = cell(N,N);
SUM_FN = zeros(N,N);
figure;
for x = 1:N,
for y = 1:N,
    plot(x,y,'ro'); hold on;
end
end
xlim([0 N+1]); ylim([0 N+1]); grid on;

for i = 1:N,
    for j = 1:N,
        P = [i,j];
        if ismember(i,[1,N]) == 1 && ismember(j,[1,N]) == 0,
            SUM_FN(i,j) = 3;
        end
        if ismember(j,[1,N]) == 1 && ismember(i,[1,N]) == 0,
            SUM_FN(i,j) = 3;
        end
        if ismember(i,[1,N]) == 1 && ismember(j,[1,N]) == 1,
            SUM_FN(i,j) = 2;
        end
        if ismember(j,[1,N]) == 0 && ismember(i,[1,N]) == 0,
            SUM_FN(i,j) = 4;
        end    
    end
end
N0 = N-1;
SUM_FN_FN = zeros(N,N);
for i = 2:N0,
    for j = 2:N0,
        P = [i,j];
        if ismember(i,[2,N0]) == 1 && ismember(j,[2,N0]) == 0,
            % ONE BOUNDARY NODE as NEIGHBOUR; All others INTERIORS
            % This Boundary node is not a corner node
            
            SUM_FN_FN(i,j) = 15;
        end
        if ismember(j,[2,N0]) == 1 && ismember(i,[2,N0]) == 0,
            SUM_FN_FN(i,j) = 15;
        end
        if ismember(i,[2,N0]) == 1 && ismember(j,[2,N0]) == 1,
            % Two BOUNDARY nodes...
            
            SUM_FN_FN(i,j) = 14;
        end
        if ismember(j,[2,N0]) == 0 && ismember(i,[2,N0]) == 0,
            % ALL INTERIOR NODES...
            
            SUM_FN_FN(i,j) = 16;
        end    
    end
end

for p = [1,N],
    for q = 1:N,
        if ismember(p,[1,N]) == 1 && ismember(q,[1,N])==1,
           SUM_FN_FN(p,q) = 6; % TWO BOUNDARY NODES
        else
            if p==2 || p == N-1,
            SUM_FN_FN(p,1) = 9; % TWO BOUNDARY NODES + ONE INTERIOR
            SUM_FN_FN(p,N) = 9;
            else
                SUM_FN_FN(p,q) = 10;
            end
            if q==2 || q == N-1,
            SUM_FN_FN(q,1) = 9; % TWO BOUNDARY NODES + ONE INTERIOR
            SUM_FN_FN(q,N) = 9;
            else
                SUM_FN_FN(q,p) = 10;
            end
             
        end
    end
end

% RANDOM WALK GENERATION
ST = C; % START NODE
% RULES for CHOOSING DESTINATION NODE
% 1) ISOLATE THE FREE NODES
% 2) DETERMINE the EXTENT of FREEDOM ASSOCIATED with these FREE NODES
% For this use the SUM_FN_FN map...
% 3) Pick the DESTINATION node for which the SUM_FN_FN SCORE is MAXIMUM

IND0 = 1:N^2;
CHAIN{1} = C;
for c = 1:N^2-1,
    Node_LOC = CHAIN{c};
    CHAIN_CODE(c) = (Node_LOC(1)) + (Node_LOC(2)-1)*N;
    
    [POT_DESTS, EMPTY] = fun_QNN(Node_LOC, N, CHAIN);
    if EMPTY == 0,
    for i = 1:length(POT_DESTS),
        A = POT_DESTS{i};
        P1 = SUM_FN(A(2), A(1));
        P2 = SUM_FN_FN(A(2), A(1));
        PP(i) = P1*P2;
    end
    [S, IND] = sort(PP,'descend');
    MAX_IND = find(PP == S(1));
    LM = length(MAX_IND);
    R = randperm(LM);      
    D = POT_DESTS{MAX_IND(R(1))};
    XL = [Node_LOC(1); D(1)];
    YL = [Node_LOC(2); D(2)];
    line(XL,YL,'Color','r','LineWidth',3,'Marker','.','LineStyle','-'); hold on;
    pause(0.5)
    CHAIN{c+1} = D;
    clear POT_DESTS PP
    else
        SET = setdiff(IND0,CHAIN_CODE);
        LS = length(SET);
        Q = randperm(LS);
        SNEW = SET(Q(1));
        Q2 = floor(SNEW/N)+1;
        Q1 = SNEW - N*(Q2-1);
        CHAIN{c+1} = [Q1,Q2];
    end
    
    %pause;
end

    
% SUM_FN
% SUM_FN_FN

% for c=1:441
%     index = CHAIN{1,c};
%     if index(1,1)==0 
%         index(1,1)= index(1,1)+1;
%     end
%     if index(1,2)==0
%         index(1,2)= index(1,2)+1;
%     end
%     encrypted_signal(1,c) = random_signal(index(1,1), index(1,2));
% end

% clear all
% close all
% I_real = imread('001_r.jpg');
% I_pp = imread('001_pp.jpg');
% I_real = mat2gray(imresize(rgb2gray(I_real),[101 101]));
% I_pp =  mat2gray(imresize(rgb2gray(I_pp),[101 101]));
% 
% [CHAIN, CHAIN_CODE, SUM_FN, SUM_FN_FN,encrypted_signal_real]= random_walk_generation(I_real,101)
% [CHAIN, CHAIN_CODE, SUM_FN, SUM_FN_FN,encrypted_signal_pp]= random_walk_generation(I_pp,101)

            