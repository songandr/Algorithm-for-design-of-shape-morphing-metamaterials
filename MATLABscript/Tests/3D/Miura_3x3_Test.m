% This script is designed as a test case for the MATLAB function 
% minimizationAlgorithm which is based on the paper:
% 
% "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 03/20/26.
%
% The initial configuration consists of 16 panels in the Miura-Ori
% configuration (4x4 unit cell). This test is interested in the folding of these 
% panels and the resultant energy calculation.

% Test Parameters
test = "axial";
crease_test = "Miura3x3";
panelRemoval = 0;

% Initial x-values
l1R = [1; 0];
l2R = [0; 1];
x1 = [0; 0];
x13 = x1 + l1R + l2R;
x4 = x1 + l1R;
x16 = x1 + l2R;
x2 = (x1+x4)/3;
x3 = 2*x2;
x8 = x16/3;
x9 = 2*x8;
x5 = x4+(x13-x4)/3;
x12 = x4+2*(x13-x4)/3;
x7 = x8+(x5-x8)/3;
x6 = x8+2*(x5-x8)/3;
x10 = x9+(x12-x9)/3;
x11 = x9+2*(x12-x9)/3;
x15 = x16+(x13-x16)/3;
x14 = x16+2*(x13-x16)/3;
gamma = pi/6;
a = -[0; sin(gamma)]/3;
x2 = x2 - a;
x7 = x7 - a;
x10 = x10 - a;
x15 = x15 - a;
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16];

% Randomized ICs
%{
rng(3, "twister");
r = normrnd(0, 0.01, [32, 1]);
% only perturb eligible DoFs
r(1:2) = 0; % fixed first node
r(end-9:end) = 0; % top nodes
r(7:10) = 0; % side nodes

x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16];
x = x+r;

% maintain original lattice vectors post-perturbation
x(9) = x(9) + r(15); % shift right side nodes
x(10) = x(10) + r(16);
x(23) = x(23) + r(17);
x(24) = x(24) + r(18);
x(end-3) = x(end-3) + r(3); % shift top side nodes
x(end-2) = x(end-2) + r(4);
x(end-5) = x(end-5) + r(5);
x(end-4) = x(end-4) + r(6);
%}
phi = zeros(length(x)/2*3, length(x));
if test == "reference"
    
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = 1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = 1;
        end
    end
    
elseif test == "axial"
    lambda_1 = 0.334; % axial deformation
    lambda_2 = 0.334;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = lambda_1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = lambda_2;
        end    
    end

elseif test == "shear"
    %{
    gamma_x = 0;

    % sheared initial x case
    phi_x = zeros(length(x), length(x));
    for i = 1:length(x)
        phi_x(i,i) = 1;
        if mod(i,2) == 1
            phi_x(i,i+1) = gamma_x;
        else
            phi_x(i,i-1) = gamma_x;
        end
    end
    x = phi_x*x;
    %}
    
    gamma = 0.20;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = 1;
            phi(i,i-floor(i/3)+1) = gamma;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = 1;
            phi(i,i-floor(i/3)-1) = gamma;
        end    
    end

elseif test == "mix"
    gamma = 0.20; % shear deformation
    lambda_1 = 0.75; % axial deformation
    lambda_2 = 0.75;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = lambda_1;
            phi(i,i-floor(i/3)+1) = gamma;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = lambda_2;
            phi(i,i-floor(i/3)-1) = gamma;
        end    
    end

elseif test == "conditioned_20"
    % xOpt from gamma = 0.20, lambda = 0.75, seed: rng(2, "twister"), 
    x = [1.26695175141596e-19; -6.48737621103359e-32; 0.267950695723525; -0.00856557126671520; 0.597275108267772; -0.0658561604695548; 0.758533148651623; -0.0477120791408657; 1; -3.18787862285775e-18; 1.06599026554045; 0.135845832613217; 0.758646412563898; 0.107940023384610; 0.640721269680140; 0.102055475820708; 0.298056417759411; 0.308396821089452; 0.0659902655404491; 0.135845832613217; -0.0627659421090146; 0.507709660968152; 0.0969766982619429; 0.648830182538741; 0.438407409494005; 0.456358508980193; 0.740652566802095; 0.486724393899420; 0.937234057890986; 0.507709660968153; 0.945282514630622; 0.807480337714779; 0.781647830773168; 0.747968088765264; 0.535941073133419; 0.722221132381985; 0.215496386781822; 0.809453296724130; -0.0547174853693774; 0.807480337714779; -5.13550312371518e-18; 1; 0.267950695723525; 0.991434428733284; 0.597275108267772; 0.934143839530445; 0.793812214934423; 0.952287920859134; 1; 1];
    % yOpt from gamma = 0.20, lambda = 0.75, seed: rng(2, "twister"), 
    y = [1.14135889874404e-33; 1.01826518992631e-18; -4.61921801435536e-18; 0.261716129567368; 0.0448276061040079; -0.0250377639217541; 0.436320144852988; 0.205026963679619; -0.259580390726072; 0.561879025144747; 0.189403475158481; -0.147345021203485; 0.75; 0.2; -1.18130607321628e-17; 0.802033348569237; 0.316962375885326; 0.0829045038932574; 0.552587446752796; 0.338116586475769; -0.0984359056810793; 0.458605529123422; 0.351880392867948; -0.171243500961822; 0.260497685677772; 0.295972764093041; 0.174172551995775; 0.0520333485692377; 0.116962375885326; 0.0829045038932573; 0.187277928653808; 0.311637804870438; -0.229870779558305; 0.338630606431490; 0.460147482377046; -0.182652037948484; 0.526791811584703; 0.544330485735345; -0.523242871832384; 0.775627234401789; 0.521549673617831; -0.346861448003352; 0.937277928653808; 0.511637804870438; -0.229870779558305; 0.917138666295388; 0.784095051681599; -0.108369649837651; 0.794389990760509; 0.757176961730687; -0.227283287020317; 0.599015100323342; 0.772288943417606; -0.374443270670963; 0.429488000194149; 0.633845970373430; -0.126630934082619; 0.167138666295388; 0.584095051681599; -0.108369649837651; 0.2; 0.75; 3.91359701774152e-19; 0.461716129567368; 0.794827606104008; -0.0250377639217542; 0.636320144852988; 0.955026963679620; -0.259580390726072; 0.788338324856848; 0.946459288415041; -0.147345021203485; 0.95; 0.95; 6.90416229466756e-18];
    
    delta = 0.07; % additional shear on top of gamma above with fixed lambdas
    for i = 1:length(y)
        if mod(i,3) == 1
            y(i) = y(i) + delta*x(i-floor(i/3)+1); % y_{i,1} = y_{0,1} + delta*x_{0,2}
            %y(i) = y(i) + delta*(x(i-floor(i/3)+1) - 0.5*x(i-floor(i/3))); % if further contraction is necessary
        elseif mod(i,3) == 2 
            y(i) = y(i) + delta*x(i-floor(i/3)-1); %y_{i,2} = y_{0,2} + delta*x_{0,1}
            %y(i) = y(i) + delta*(x(i-floor(i/3)-1) - 0.5*x(i-floor(i/3))); % if further contraction is necessary
        end
    end

elseif test == "conditioned_28"
    % xOpt from gamma = 0.28, lambda = 0.75; see conditioned_20 above 
    x = [1.16956000141296e-19; -6.85959494041963e-32; 0.263648456054710; -0.00830938488824884; 0.638133105865052; -0.0878169176446899; 0.773635949695853; -0.0662779093887950; 1; -3.34656760387148e-18; 1.07923943639164; 0.125108796340599; 0.774551511425953; 0.0625961887884249; 0.688748415828559; 0.0456780964248715; 0.300553839922697; 0.321186544895071; 0.0792394363916433; 0.125108796340599; -0.106776068316943; 0.546769856682162; 0.0429859317595547; 0.701721561987143; 0.429125786001006; 0.439880277027983; 0.733143035915184; 0.507421766998406; 0.893223931683058; 0.546769856682163; 0.929141910635002; 0.823871914370499; 0.782699927534252; 0.748967269672129; 0.549151002416303; 0.717073009587014; 0.177110013269285; 0.837588870875085; -0.0708580893649972; 0.823871914370499; -7.04052434619803e-18; 1; 0.263648456054710; 0.991690615111751; 0.638133105865052; 0.912183082355310; 0.808915015978653; 0.933722090611205; 1; 1];
    % yOpt from gamma = 0.28, lambda = 0.75; see conditioned_20 above 
    y = [-3.48915169333884e-33; 1.77392444356641e-18; -7.58024365339307e-18; 0.251017434769753; 0.0694618453458889; -0.0321646715760294; 0.485813659054501; 0.265056611053120; -0.261151469097747; 0.585189905934480; 0.252205277795856; -0.156610439161887; 0.750000000000000; 0.280000000000000; -1.23328617921278e-17; 0.812095854224971; 0.385395798070430; 0.0865784401108210; 0.581426721746624; 0.378130651691637; -0.121459633912567; 0.515361219953316; 0.378995872252953; -0.180968159094347; 0.246597245576512; 0.301866974286361; 0.207605813524479; 0.0620958542249716; 0.105395798070430; 0.0865784401108209; 0.239920942453258; 0.359680540331123; -0.252544301276526; 0.371247268463906; 0.520175506003352; -0.183191798458047; 0.626221864416010; 0.626352405104659; -0.566135920182750; 0.864150887688107; 0.630410266180165; -0.361916776453368; 0.989920942453258; 0.639680540331123; -0.252544301276526; 0.985621390514721; 0.881098446679265; -0.111516297288067; 0.882459189039047; 0.840760102588080; -0.233272844606698; 0.703335975439426; 0.852445529748230; -0.383765372698334; 0.470143549194844; 0.677304705828165; -0.125419218810888; 0.235621390514721; 0.601098446679265; -0.111516297288067; 0.280000000000000; 0.750000000000000; 8.91667578561599e-19; 0.531017434769753; 0.819461845345889; -0.0321646715760295; 0.765813659054501; 1.01505661105312; -0.261151469097747; 0.891649205646581; 1.01208341635504; -0.156610439161887; 1.03000000000000; 1.03000000000000; 9.89984209184458e-18];
    
    delta = 0.07; % additional shear on top of gamma above with fixed lambdas
    for i = 1:length(y)
        if mod(i,3) == 1
            %y(i) = y(i) + delta*x(i-floor(i/3)+1); % y_{i,1} = y_{0,1} + delta*x_{0,2}
            y(i) = y(i) + delta*(x(i-floor(i/3)+1) - 1.1*x(i-floor(i/3))); % if further contraction is necessary
        elseif mod(i,3) == 2 
            %y(i) = y(i) + delta*x(i-floor(i/3)-1); %y_{i,2} = y_{0,2} + delta*x_{0,1}
            y(i) = y(i) + delta*(x(i-floor(i/3)-1) - 1.1*x(i-floor(i/3))); % if further contraction is necessary
        end
    end

elseif test == "conditioned_35"
    % xOpt from gamma = 0.35, lambda = 0.673; see conditioned_28 above 
    x = [7.07759081486666e-20; -5.51893371254605e-32; 0.254334645538455; -0.00452203627420662; 0.642724914046687; -0.0909264406223934; 0.772487897978393; -0.0679042150824031; 1; -2.65929576639134e-18; 1.07618058031525; 0.125305403182899; 0.771715759718018; 0.0524912480065814; 0.691171801446607; 0.0345882632002104; 0.298602394889495; 0.323147343770554; 0.0761805803152487; 0.125305403182899; -0.132735150538694; 0.564495636985668; 0.0231261846951282; 0.718323801519067; 0.419778688222880; 0.442236006210470; 0.721583439188999; 0.522102778860930; 0.867264849461307; 0.564495636985669; 0.922354384584442; 0.831285987513876; 0.777376625639663; 0.753257773533795; 0.546285424557036; 0.719781727197889; 0.157832505162507; 0.851517996840257; -0.0776456154155572; 0.831285987513876; -6.63254208647857e-18; 1; 0.254334645538455; 0.995477963725793; 0.642724914046687; 0.909073559377607; 0.807766964261193; 0.932095784917597; 1; 1];
    % yOpt from gamma = 0.35, lambda = 0.673; see conditioned_28 above 
    y = [-9.00561201088777e-21; 2.13058468197095e-18; -9.06247990314967e-18; 0.228630730354827; 0.0993689910724733; -0.0437588550454732; 0.433472436697156; 0.344256134628057; -0.280360512091404; 0.521185351527859; 0.326188632730465; -0.170228533031100; 0.673000000000000; 0.350000000000000; -1.37918750133468e-17; 0.728493561643764; 0.453178507949612; 0.0914440710247845; 0.514687138493365; 0.443434542799739; -0.135728097554631; 0.458242433780564; 0.445396460297171; -0.196603527568431; 0.229237004035812; 0.310864666731938; 0.213539416349996; 0.0554935616437648; 0.103178507949613; 0.0914440710247844; 0.297070691335527; 0.311902963831328; -0.274198975353037; 0.422449726460631; 0.478760304060805; -0.199077219498326; 0.639431327641820; 0.644894136285250; -0.604194592832128; 0.861836989418320; 0.651615423935102; -0.382431984286410; 0.970070691335527; 0.661902963831328; -0.274198975353037; 0.974234087625518; 0.885152165303301; -0.117265033069350; 0.879968735861096; 0.846872712641261; -0.247102676288572; 0.714029309471926; 0.864488380964775; -0.408390551584954; 0.511816488804132; 0.637070971153277; -0.135997735022618; 0.301234087625518; 0.535152165303301; -0.117265033069350; 0.350000000000000; 0.673000000000000; 1.02378131830864e-18; 0.578630730354827; 0.772368991072473; -0.0437588550454733; 0.783472436697156; 1.01725613462806; -0.280360512091404; 0.894928163136184; 1.01153630592944; -0.170228533031100; 1.02300000000000; 1.02300000000000; 1.17115599686443e-17];

    delta = 0.05; % additional shear on top of gamma above with fixed lambdas
    for i = 1:length(y)
        if mod(i,3) == 1
            %y(i) = y(i) + delta*x(i-floor(i/3)+1); % y_{i,1} = y_{0,1} + delta*x_{0,2}
            y(i) = y(i) + delta*(x(i-floor(i/3)+1) - 0.5*x(i-floor(i/3))); % if further contraction is necessary
        elseif mod(i,3) == 2 
            %y(i) = y(i) + delta*x(i-floor(i/3)-1); %y_{i,2} = y_{0,2} + delta*x_{0,1}
            y(i) = y(i) + delta*(x(i-floor(i/3)-1) - 0.5*x(i-floor(i/3))); % if further contraction is necessary
        end
    end

end

if contains(test, "conditioned") == 0 % only applicable for ICs at the original reference config
    y = phi*x;
    % Bias mountain y's in the z-axis
    epsilon = 0.05;
    y(15) = y(15)+epsilon;
    y(18) = y(18)+epsilon;
    y(21) = y(21)+epsilon;
    y(24) = y(24)+epsilon;
    %{
    y(39) = y(39)+epsilon;
    y(42) = y(42)+epsilon;
    y(45) = y(45)+epsilon;
    y(48) = y(48)+epsilon;
    %}
end

% x rigidity constraint matrix (9x16) of (2x2)
L_0 = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), -eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2)];
    
% y rigidity constraint matrix (9x16) of (3x3)
L = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), -eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3)];

% populate the vector numbering all of the panels
J = 1:9;

% index set for each panel 
F1 = [1, 2, 7, 8];
F2 = [2, 3, 6, 7];
F3 = [3, 4, 5, 6];
F4 = [5, 6, 11, 12];
F5 = [6, 7, 10, 11];
F6 = [7, 8, 9, 10];
F7 = [9, 10, 15, 16];
F8 = [10, 11, 14, 15];
F9 = [11, 12, 13, 14];

T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;
T5 = F5;
T6 = F6;
T7 = F7;
T8 = F8;
T9 = F9;

% 3D array containing index set of x coordinates for panel j
Tj = cell(length(J), 1);
for j=1:length(J)
    Tj{j} = eval(sprintf('T%d', j));
end

Fj = Tj;

% Initial R 
R = cell(length(J), 1);
for j = 1:length(J)
    R{j} = eye(3); % identity matrix
end

% initial tolerance for minimization
tol = 10^(-5);

[yOpt, xOpt, Ropt] = minimizationAlgorithm(x, y, Fj, Tj, J, R, L, L_0, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;
plot4vectors3D(vectors, titles, visualizeLatticeVec, crease_test);

%{
% post-processing
% random perturbation post-processing step
E = 1; % initialize energy
postprocessing_count = 0;
while E > 5*10^(-5)
r = normrnd(0, 0.001, [length(x), 1]);
r_y = normrnd(0, 0.001, [length(y), 1]);

% only perturb eligible DoFs
r(1:2) = 0; % fixed first node
r(41:50) = 0; % top nodes
r(9:12) = 0; % side nodes
r(29:32) = 0; % side nodes
r_y(1:3) = 0; % fixed first node
r_y(end-14:end) = 0; % top nodes
r_y(13:18) = 0; % side nodes
r_y(43:48) = 0; % side nodes

x_perturbed = xOpt + r;
y_perturbed = yOpt + r_y;

% maintain original lattice vectors post-perturbation
x_perturbed(11) = xOpt(11) + r(19); % shift right side nodes
x_perturbed(12) = xOpt(12) + r(20);
x_perturbed(29) = xOpt(29) + r(21);
x_perturbed(30) = xOpt(30) + r(22);
x_perturbed(31) = xOpt(31) + r(39);
x_perturbed(32) = xOpt(32) + r(40);
x_perturbed(43) = xOpt(43) + r(3); % shift top side nodes
x_perturbed(44) = xOpt(44) + r(4);
x_perturbed(45) = xOpt(45) + r(5);
x_perturbed(46) = xOpt(46) + r(6);
x_perturbed(47) = xOpt(47) + r(7);
x_perturbed(48) = xOpt(48) + r(8);

y_perturbed(16) = yOpt(16) + r_y(28); % shift right side nodes
y_perturbed(17) = yOpt(17) + r_y(29);
y_perturbed(18) = yOpt(18) + r_y(30);
y_perturbed(43) = yOpt(43) + r_y(31);
y_perturbed(44) = yOpt(44) + r_y(32);
y_perturbed(45) = yOpt(45) + r_y(33);
y_perturbed(46) = yOpt(46) + r_y(58);
y_perturbed(47) = yOpt(47) + r_y(59);
y_perturbed(48) = yOpt(48) + r_y(60);
y_perturbed(64) = yOpt(64) + r_y(4); % shift top side nodes
y_perturbed(65) = yOpt(65) + r_y(5);
y_perturbed(66) = yOpt(66) + r_y(6);
y_perturbed(67) = yOpt(67) + r_y(7);
y_perturbed(68) = yOpt(68) + r_y(8);
y_perturbed(69) = yOpt(69) + r_y(9);
y_perturbed(70) = yOpt(70) + r_y(10);
y_perturbed(71) = yOpt(71) + r_y(11);
y_perturbed(72) = yOpt(72) + r_y(12);

disp("Starting perturbation minimization...")
[yOpt, xOpt, Ropt] = minimizationAlgorithm(x_perturbed, y_perturbed, Fj, Tj, J, Ropt, L, L_0, tol);

% Compute energy associated with solution
E = 0;

% Construct necessary cj and rij vectors
cj = cell(length(J));
rij = cell(length(J));
for j = 1:length(J)
    % center of the panel calculation based on initial y vector
    [cj{j}, ~] = centerOfPanel3D(Fj{j}, yOpt);

    % pos vectors with respect to the center of the panel
    [~, rij{j}] = centerOfPanel2D(Tj{j}, xOpt);
end

for j = 1:length(J)
    for i = 1:length(Fj{j})
        k = Fj{j}(i);

        rij_temp = rij{j}(2*i-1:2*i);
        rij1 = rij_temp(1);
        rij2 = rij_temp(2);
        
        E = E + norm(yOpt(3*k-2:3*k, 1) - cj{j} - Ropt{j}*[rij1; rij2; 0])^2;
    end
end 
postprocessing_count = postprocessing_count + 1;
end

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;
plot4vectors3D(vectors, titles, visualizeLatticeVec, crease_test);

% check for concavity
if (xOpt(24)-xOpt(22))/(xOpt(23)-xOpt(21)) * (xOpt(37)-xOpt(21)) + xOpt(22) < xOpt(38)
    disp("Convex!")
else
    disp("Concave...")
end
disp("Number of post-processing steps required: "+postprocessing_count)
%}