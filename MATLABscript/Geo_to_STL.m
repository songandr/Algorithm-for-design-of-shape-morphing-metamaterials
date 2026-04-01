close all; 
clear all; 
clc;
%%

Node_ref = [3.42E-19
-1.77E-31
0.218343074
0.122522272
0.542145667
0.148918339
0.756551751
-0.081590456
1
-8.71E-18
0.998894389
0.217359992
0.721604314
0.126925908
0.506571537
0.35762545
0.200418635
0.332016095
-0.001105611
0.217359992
-0.023714142
0.538754176
0.244907297
0.509356373
0.498410398
0.557174406
0.742235752
0.470042327
0.976285858
0.538754176
0.981457844
0.776706913
0.753254663
0.763615208
0.499493594
0.791598719
0.287012314
0.806219107
-0.018542156
0.776706913
5.42E-18
1
0.218343074
1.122522272
0.542145667
1.148918339
0.756551751
0.918409544
1
1];

Node_def=[1.05E-33
1.36E-18
-5.04E-18
0.200729817
0.137070228
0.054711214
0.522109971
0.167117342
0.070291537
0.663993313
-0.090783126
-0.038415677
0.9
2.88E-17
-6.04E-18
0.894079387
0.135893924
0.170176906
0.625833138
0.033960958
0.128478421
0.483204974
0.29236519
0.236954321
0.179272046
0.263786682
0.221121891
-0.005920613
0.135893924
0.170176906
-0.024624954
0.447244258
0.088625743
0.22329625
0.440064871
0.196930104
0.476975263
0.489534657
0.205431648
0.646941651
0.369285657
0.050898931
0.875375046
0.447244258
0.088625743
0.882286858
0.679713995
0.037796614
0.65887185
0.657052298
-0.009992491
0.479254731
0.720606324
0.161623583
0.26614377
0.732496884
0.144304074
-0.017713142
0.679713995
0.037796614
4.20E-18
0.9
-1.43E-18
0.200729817
1.037070228
0.054711214
0.522109971
1.067117342
0.070291537
0.663993313
0.809216874
-0.038415677
0.9
0.9
8.71E-18

    ];

%% Orth to Non-Orth
% Node_ref = xOpt;Node_def = yOpt;

%% 
Node_ref = reshape(Node_ref, 2, [])';
Node_def = reshape(Node_def, 3, [])';

Connectivity = [
    1 2 9 10;
    2 3 8 9;
    3 4 7 8;
    4 5 6 7;
    10 9 12 11;
    9 8 13 12;
    8 7 14 13;
    7 6 15 14;
    11 12 19 20;
    12 13 18 19;
    13 14 17 18;
    14 15 16 17;
    20 19 22 21;
    19 18 23 22;
    18 17 24 23;
    17 16 25 24
    ];

figure('Color','w');
tiledlayout(1,2);

%% Reference configuration
nexttile;

Node_ref_3D = [Node_ref, zeros(size(Node_ref,1),1)];

patch('Vertices', Node_ref_3D, ...
      'Faces', Connectivity, ...
      'FaceColor', [0.8 0.9 1.0], ...
      'EdgeColor', 'k', ...
      'FaceAlpha', 0.9);

axis equal;
grid on;
view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Reference Configuration');
camlight headlight;
lighting gouraud;
hold on;

for i = 1:size(Node_ref_3D,1)
    text(Node_ref_3D(i,1), Node_ref_3D(i,2), Node_ref_3D(i,3), ...
        sprintf('  %d', i), 'FontSize', 8, 'Color', 'r');
end

%% Deformed configuration
nexttile;

patch('Vertices', Node_def, ...
      'Faces', Connectivity, ...
      'FaceColor', [0.8 0.9 1.0], ...
      'EdgeColor', 'k', ...
      'FaceAlpha', 0.9);

axis equal;
grid on;
view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Deformed Configuration');
camlight headlight;
lighting gouraud;
hold on;

for i = 1:size(Node_def,1)
    text(Node_def(i,1), Node_def(i,2), Node_def(i,3), ...
        sprintf('  %d', i), 'FontSize', 8, 'Color', 'r');
end

%% Export to STL
triFaces = [
    Connectivity(:,[1 2 3]);
    Connectivity(:,[1 3 4])
];

stlFileName_ref = 'miura_design_ref.stl';
Node_ref_3D = [Node_ref, zeros(size(Node_ref,1),1)];
write_ascii_stl(stlFileName_ref, triFaces, Node_ref_3D);
disp(['Reference STL file written to: ', stlFileName_ref]);
stlFileName_def = 'miura_design_def.stl';
write_ascii_stl(stlFileName_def, triFaces, Node_def);
disp(['STL file written to: ', stlFileName_def]);



%% Helper function
function write_ascii_stl(filename, faces, vertices)
    fid = fopen(filename, 'w');
    if fid < 0
        error('Could not open STL file for writing.');
    end

    fprintf(fid, 'solid miura_design\n');

    for f = 1:size(faces,1)
        v1 = vertices(faces(f,1),:);
        v2 = vertices(faces(f,2),:);
        v3 = vertices(faces(f,3),:);

        n = cross(v2 - v1, v3 - v1);
        if norm(n) > 0
            n = n / norm(n);
        else
            n = [0 0 0];
        end

        fprintf(fid, '  facet normal %.6e %.6e %.6e\n', n(1), n(2), n(3));
        fprintf(fid, '    outer loop\n');
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', v1(1), v1(2), v1(3));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', v2(1), v2(2), v2(3));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', v3(1), v3(2), v3(3));
        fprintf(fid, '    endloop\n');
        fprintf(fid, '  endfacet\n');
    end

    fprintf(fid, 'endsolid miura_design\n');
    fclose(fid);
end