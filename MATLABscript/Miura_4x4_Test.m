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
% Updated Date: 06/12/26.
%
% The initial configuration consists of 16 panels in the Miura-Ori
% configuration (4x4 unit cell). This test is interested in the folding of these 
% panels and the resultant energy calculation.
close all;
clear all; 
clc;

% Test Parameters
test = "mix";
crease_test = "Miura4x4"; % standard 4x4 grid crease pattern
%crease_test = "Miura4x4_diagonal1"; % for alternative crease patterns
panelRemoval = 1; % must be 1 or 2 for shear/mix tests
check_strain = 0;
slit_test = 0;

% Initial x-values
l1R = [1; 0];
l2R = [0; 1];
x1 = [0; 0];
x5 = x1 + l1R;
x21 = x1 + l2R;
x25 = x1 + l1R + l2R;
x3 = (x1+x5)/2;
x11 = (x1+x21)/2;
x15 = (x5+x25)/2;
x13 = (x11+x15)/2;
x23 = (x21+x25)/2;
x2 = (x1+x3)/2;
x4 = (x3+x5)/2;
x6 = (x5+x15)/2;
x8 = (x3+x13)/2;
x7 = (x6+x8)/2;
x10 = (x1+x11)/2;
x9 = (x8+x10)/2;
x12 = (x11+x13)/2;
x14 = (x13+x15)/2;
x16 = (x15+x25)/2;
x18 = (x13+x23)/2;
x17 = (x16+x18)/2;
x20 = (x11+x21)/2;
x19 = (x18+x20)/2;
x22 = (x21+x23)/2;
x24 = (x23+x25)/2;

% test
%{
delta = [0; 0.07];
x2 = x3 - l1R/4 + delta;
x4 = x5 - l1R/4 - delta;
x10 = (3*x1 + 2*x11)/5;
x9 = x10 + l1R/4 + delta;
x12 = x11 + l1R/4 - delta;
x14 = x13 + l1R/4 - delta;
x8 = x10 + l1R/2;
x7 = x8 + l1R/4 + delta;
x6 = x10 + l1R;
x22 = x2 + l2R;
x23 = x3 + l2R;
x24 = x4 + l2R;
%}
% Randomized ICs
rng(6, "twister");
r = normrnd(0, 1/50, [50, 1]);
% only perturb eligible DoFs
r(1:2) = 0; % fixed first node
r(41:50) = 0; % top nodes
r(9:12) = 0; % side nodes
r(29:32) = 0; % side nodes

x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25];

x = x+r;

% maintain original lattice vectors post-perturbation
x(11) = x(11) + r(19); % shift right side nodes
x(12) = x(12) + r(20);
x(29) = x(29) + r(21);
x(30) = x(30) + r(22);
x(31) = x(31) + r(39);
x(32) = x(32) + r(40);
x(43) = x(43) + r(3); % shift top side nodes
x(44) = x(44) + r(4);
x(45) = x(45) + r(5);
x(46) = x(46) + r(6);
x(47) = x(47) + r(7);
x(48) = x(48) + r(8);

% slit test: initialize fixed slit shape
if slit_test
    % set nodes 8, 9, 12 using:
    % random perturbations
    random = 1;
    
    if ~random % known solution for gamma=0.2, lambda=0.75
        x(15:16) = [0.650888890086792; 0.104143712602226];
        x(17:18) = [0.279257640210828; 0.299430803195907];
        x(23:24) = [0.100483120971925; 0.670613637764285];
    end
    
    x(25:26) = x(15:16) + (x(23:24) - x(17:18)); % make parallelogram
    if panelRemoval == 2
        x(33:34) = x(27:28) + (x(35:36) - x(25:26)); % make 2nd parallelogram
    end
end

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
    lambda_1 = 0.8; % axial deformation
    lambda_2 = 0.6;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = lambda_1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = lambda_2;
        end    
    end

elseif test == "axial_conditioned_0.9"
    % xOpt, yOpt from lambda_1 = lambda_2 = 0.9, seed: rng(3, "twister"), 20 post processing perturbations 
    x = [3.41518115915436e-19;-1.77434731347731e-31;0.218343073588351;0.122522272106046;0.542145666615869;0.148918338874566;0.756551751364812;-0.0815904562051213;1;-8.71196966144387e-18;0.998894388778316;0.217359991897539;0.721604314204667;0.126925907858072;0.506571536504127;0.357625450095033;0.200418635143094;0.332016094813508;-0.00110561122168521;0.217359991897539;-0.0237141420178077;0.538754176016250;0.244907297013311;0.509356373280514;0.498410397749940;0.557174405790517;0.742235752030443;0.470042326789816;0.976285857982193;0.538754176016251;0.981457844265286;0.776706912982461;0.753254663041508;0.763615207612582;0.499493593936985;0.791598718679399;0.287012314311251;0.806219107287940;-0.0185421557347141;0.776706912982461;5.42216028917731e-18;1;0.218343073588351;1.12252227210605;0.542145666615869;1.14891833887457;0.756551751364812;0.918409543794878;1;1];
    y = [1.05051775048501e-33;1.36480613528626e-18;-5.03660726210027e-18;0.200729816582132;0.137070227815271;0.0547112142286052;0.522109970957241;0.167117341731585;0.0702915370704601;0.663993313470627;-0.0907831263018998;-0.0384156766683944;0.900000000000000;2.87882107635592e-17;-6.03981294511484e-18;0.894079386522476;0.135893924421486;0.170176906133031;0.625833138159709;0.0339609578671692;0.128478420869026;0.483204974429579;0.292365189526130;0.236954321083634;0.179272045726360;0.263786681965989;0.221121891072435;-0.00592061347752497;0.135893924421486;0.170176906133031;-0.0246249537983399;0.447244257724843;0.0886257432363102;0.223296249892288;0.440064870943731;0.196930104098718;0.476975263376644;0.489534656866602;0.205431648041772;0.646941650936278;0.369285656600793;0.0508989314778936;0.875375046201660;0.447244257724843;0.0886257432363102;0.882286857670766;0.679713994724449;0.0377966138896006;0.658871849872141;0.657052298250449;-0.00999249142339743;0.479254730753498;0.720606323874594;0.161623583217648;0.266143770031365;0.732496883574956;0.144304073605431;-0.0177131423292340;0.679713994724449;0.0377966138896005;4.19804585172537e-18;0.900000000000000;-1.43430097110850e-18;0.200729816582132;1.03707022781527;0.0547112142286053;0.522109970957241;1.06711734173158;0.0702915370704601;0.663993313470627;0.809216873698101;-0.0384156766683943;0.900000000000000;0.900000000000000;8.70675169969928e-18];
    % additional axial deformation on y
    lambda_1 = 0.8;
    lambda_2 = 0.8;
    
    phi = zeros(length(y));
    for i = 1:length(y)
        if mod(i,3) == 1
            phi(i,i) = lambda_1/0.9;
        elseif mod(i,3) == 2
            phi(i,i) = lambda_2/0.9;
        else
            phi(i,i) = 1;
        end    
    end
    y = phi*y;
    
elseif test == "shear"
    % sheared initial x case
    %{
    gamma_x = 0;
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
    gamma = 0.2; % shear deformation
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
    % new xOpt using post processing clean up methods
    x = [-5.63251315268542e-19;1.67820410905161e-31;0.139494982385109;-0.103702143084424;0.521732036052458;-0.216136702675179;0.738780186820239;-0.246735839673419;1;8.58198579747825e-18;0.974195844083898;0.0478881611374225;0.605764099656293;-0.117648706909661;0.559527674045171;-0.141869173036313;0.319313697704225;0.235351596946954;-0.0258041559161092;0.0478881611374222;-0.122274254694976;0.561052545516764;0.0664446293894950;0.660385048000903;0.320250007940977;0.295566027214427;0.641049896008555;0.448620668104801;0.877725745305024;0.561052545516765;0.967667187886938;0.753930663845213;0.735229385007454;0.695543172233824;0.469071735176112;0.601847707527864;0.0669335582174872;0.661692151107075;-0.0323328121130586;0.753930663845213;5.19405640878180e-18;1;0.139494982385109;0.896297856915576;0.521732036052458;0.783863297324816;0.774059253103038;0.753264160326583;1;1];
    % yOpt from gamma = 0.28, lambda = 0.75; see conditioned_20 above 
    y = [-3.48915169333884e-33; 1.77392444356641e-18; -7.58024365339307e-18; 0.251017434769753; 0.0694618453458889; -0.0321646715760294; 0.485813659054501; 0.265056611053120; -0.261151469097747; 0.585189905934480; 0.252205277795856; -0.156610439161887; 0.750000000000000; 0.280000000000000; -1.23328617921278e-17; 0.812095854224971; 0.385395798070430; 0.0865784401108210; 0.581426721746624; 0.378130651691637; -0.121459633912567; 0.515361219953316; 0.378995872252953; -0.180968159094347; 0.246597245576512; 0.301866974286361; 0.207605813524479; 0.0620958542249716; 0.105395798070430; 0.0865784401108209; 0.239920942453258; 0.359680540331123; -0.252544301276526; 0.371247268463906; 0.520175506003352; -0.183191798458047; 0.626221864416010; 0.626352405104659; -0.566135920182750; 0.864150887688107; 0.630410266180165; -0.361916776453368; 0.989920942453258; 0.639680540331123; -0.252544301276526; 0.985621390514721; 0.881098446679265; -0.111516297288067; 0.882459189039047; 0.840760102588080; -0.233272844606698; 0.703335975439426; 0.852445529748230; -0.383765372698334; 0.470143549194844; 0.677304705828165; -0.125419218810888; 0.235621390514721; 0.601098446679265; -0.111516297288067; 0.280000000000000; 0.750000000000000; 8.91667578561599e-19; 0.531017434769753; 0.819461845345889; -0.0321646715760295; 0.765813659054501; 1.01505661105312; -0.261151469097747; 0.891649205646581; 1.01208341635504; -0.156610439161887; 1.03000000000000; 1.03000000000000; 9.89984209184458e-18];
    % new yOpt using post processing clean up methods
    y = [-1.10681393457758e-33;4.30325965056969e-18;-1.75727022332155e-17;0.141526424814643;-0.100496777044577;-0.00856760034187486;0.328495991743004;0.125462677132323;-0.277937883417519;0.542310950820319;0.0827573075925383;-0.222350807351853;0.750000000000000;0.270000000000000;-7.27353685527883e-18;0.721302979245862;0.313758971350425;0.0160788587004369;0.400868203568531;0.203114859078997;-0.203274342159298;0.359088497745558;0.186162541712691;-0.229199673054434;0.283538557637639;0.179782510108695;0.211910548329753;-0.0286970207541388;0.0437589713504213;0.0160788587004367;0.208265354337544;0.399225794255763;-0.284605819501074;0.376832412758749;0.471594249378657;-0.176213063182225;0.467022239148736;0.488350103054510;-0.610982622538714;0.749502979418140;0.593279904247951;-0.422815573822111;0.958265354337544;0.669225794255763;-0.284605819501075;1.02327423781014;0.828436721029936;-0.158909790151300;0.812481855471685;0.799559205805490;-0.269567120404470;0.574202192904441;0.740987957466006;-0.409180442653646;0.377133170844361;0.472764771365377;-0.175438625736366;0.273274237810136;0.558436721029935;-0.158909790151300;0.270000000000000;0.750000000000000;3.82814317032276e-18;0.411526424814642;0.649503222955418;-0.00856760034187490;0.598495991743000;0.875462677132330;-0.277937883417519;0.838770250532421;0.842282655488899;-0.222350807351853;1.02000000000000;1.02000000000000;1.96058271741934e-17];
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
    y(18) = y(18)+epsilon;
    y(21) = y(21)+epsilon;
    y(24) = y(24)+epsilon;
    y(27) = y(27)+epsilon;
    y(30) = y(30)+epsilon;
    y(48) = y(48)+epsilon;
    y(51) = y(51)+epsilon;
    y(54) = y(54)+epsilon;
    y(57) = y(57)+epsilon;
    y(60) = y(60)+epsilon;
end

% x rigidity constraint matrix (10x25) of (2x2)
L_0 = [-eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
    
% y rigidity constraint matrix (10x25) of (3x3)
L = [-eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3);
    -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)];

% fixed slit constraint for slit test
if slit_test
if panelRemoval==1 % one parallelogram slit
L_0 = [zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), -eye(2), zeros(2), zeros(2), eye(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
end
if panelRemoval==2 % two parallelogram slits
L_0 = [zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), -eye(2), zeros(2), zeros(2), eye(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), eye(2), zeros(2), zeros(2), -eye(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);
     eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
end
end

% populate the vector numbering all of the panels
J = 1:16;

% index set for each panel 
F1 = [1, 2, 9, 10];
F2 = [2, 3, 8, 9];
F3 = [3, 4, 7, 8];
F4 = [4, 5, 6, 7];
F5 = [6, 7, 14, 15];
F6 = [7, 8, 13, 14];
F7 = [8, 9, 12, 13];
F8 = [9, 10, 11, 12];
F9 = [11, 12, 19, 20];
F10 = [12, 13, 18, 19];
F11 = [13, 14, 17, 18];
F12 = [14, 15, 16, 17];
F13 = [16, 17, 24, 25];
F14 = [17, 18, 23, 24];
F15 = [18, 19, 22, 23];
F16 = [19, 20, 21, 22];

if panelRemoval == 1
    J = 1:15;
    F7 = F8;
    F8 = F9;
    F9 = F10;
    F10 = F11;
    F11 = F12;
    F12 = F13;
    F13 = F14;
    F14 = F15;
    F15 = F16;
elseif panelRemoval == 2
    J = 1:14;
    F7 = F8;
    F8 = F9;
    F9 = F10;
    F10 = F12;
    F11 = F13;
    F12 = F14;
    F13 = F15;
    F14 = F16;
end

T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;
T5 = F5;
T6 = F6;
T7 = F7;
T8 = F8;
T9 = F9;
T10 = F10;
T11 = F11;
T12 = F12;
T13 = F13;
T14 = F14;
if panelRemoval ~= 2
    T15 = F15;
end
if panelRemoval == 0
    T16 = F16;
end

% 3D array containing index set of x coordinates for panel j
Tj = cell(length(J), 1);
for j=1:length(J)
    Tj{j} = eval(sprintf('T%d', j));
end

Fj = Tj;

% Initial R 
for j = 1:length(J)
    R{j} = eye(3); % identity matrix
end

% initial tolerance for minimization
tol = 10^(-5);

figure
ax = gca;
hold(ax,'on')
set(ax,'YScale','log')
ax.FontSize = 16;
[yOpt, xOpt, Ropt, ax, iterOffset] = minimizationAlgorithm(x, y, Fj, Tj, J, R, L, L_0, tol, ax);

% post-processing
E = 1; % initialize energy
count = 0;
max_count = 1000;
while E > tol && count < max_count
    % random perturbation post-processing step
    rng(5, "twister")
    r = normrnd(0, 10^(-5), [length(x), 1]);
    r_y = normrnd(0, 10^(-5), [length(y), 1]);
    
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
    if slit_test
    x_perturbed(25:26) = xOpt(25:26) + r(15:16) - r(17:18) + r(23:24); % shift x13 to enforce parallelogram slit
    if panelRemoval == 2
    x_perturbed(33:34) = xOpt(33:34) + r(27:28) - r(25:26) + r(35:36); % shift x17 to enforce parallelogram slit
    end
    end
    
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
    [yOpt_perturbed, xOpt_perturbed, Ropt, ax, iterOffset] = minimizationAlgorithm(x_perturbed, y_perturbed, Fj, Tj, J, R, L, L_0, tol, ax, iterOffset);

    xOpt = xOpt_perturbed;
    yOpt = yOpt_perturbed;

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
    disp("Energy: " + E)
    if E < tol || count == max_count-1
        titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
        vectors = {x, y, xOpt, yOpt};
        visualizeLatticeVec = true;
        plot4vectors3D(vectors, titles, visualizeLatticeVec, crease_test);
    end
    count = count + 1;
end

disp("Total post-processing steps: "+count)

% check for concavity
if (xOpt(24)-xOpt(22))/(xOpt(23)-xOpt(21)) * (xOpt(37)-xOpt(21)) + xOpt(22) < xOpt(38)
    disp("Convex!")
else
    disp("Concave...")
end

% check strains
if check_strain
    nodes = reshape(1:25, 5, 5)'; % Standard row-major grid
    nodes(2:2:end, :) = fliplr(nodes(2:2:end, :)); % Flip even rows for snake pattern
    E_h = [reshape(nodes(:, 1:end-1), [], 1), reshape(nodes(:, 2:end), [], 1)]; % (connect neighbors in the same row)
    E_v = [reshape(nodes(1:end-1, :), [], 1), reshape(nodes(2:end, :), [], 1)]; % (connect neighbors in the same column)
    edges = [E_h; E_v]; % combine
    edges = sort(edges, 2); % [smaller_node, larger_node] per row
    edges = sortrows(unique(edges, 'rows')); % remove duplicates, sort the list
    l0 = zeros(length(edges), 1);
    l = zeros(length(edges), 1);
    strains = zeros(length(edges), 1);
    % compute strains for each edge
    for i = 1:length(edges)
        x_a = [xOpt(2*edges(i,1)-1), xOpt(2*edges(i,1))]; % node 1 of the i-th edge
        x_b = [xOpt(2*edges(i,2)-1), xOpt(2*edges(i,2))]; % node 2 of the i-th edge
        l0(i) = norm(x_a - x_b); % reference configuration lengths
        y_a = [yOpt(3*edges(i,1)-2), yOpt(3*edges(i,1)-1), yOpt(3*edges(i,1))];
        y_b = [yOpt(3*edges(i,2)-2), yOpt(3*edges(i,2)-1), yOpt(3*edges(i,2))];
        l(i) = norm(y_a - y_b); % deformed configuration lengths
        strains(i) = (l(i) - l0(i))/l0(i); % resultant strains from the deformation
    end
    [max_val, max_idx] = max(abs(strains));
    disp("Max strain: " + strains(max_idx) + " @ edge (" + edges(max_idx, 1) + ", " + edges(max_idx, 2) + ")")
    disp("Average absolute strain: " + mean(abs(strains)))
    disp("Total absolute strain: " + sum(abs(strains)))
end


% compute geometric stiffness
K = geometricStiffness(xOpt, yOpt, Ropt, Tj, Fj, J, L, 5, 4, 3, "translation");
