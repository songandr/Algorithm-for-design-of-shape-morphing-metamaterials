function plot1vector(vec, name, title1, vis, test)

% plots a column vector of 3-dimensional points (x, y, z1, x2, y2,
% z2....)

% vis is a boolean that determines whether the visualization of the
% symmetry vectors is off or on

gap = 0.4;
vec = cell2mat(vec);
x = zeros(length(vec)/2, 1);
y = zeros(length(vec)/2, 1);

for i = 1:length(vec)/2
x(i) = vec(2*i-1);
y(i) = vec(2*i);
labels{i} = [name, num2str(i)];
end

figure
title(title1)
xlabel('x');
ylabel('y');

hold on
axis equal
xlim([min(x) - gap, max(x) + gap]);
ylim([min(y) - gap, max(y) + gap]);
plot(x, y);

for i = 1:length(labels)
    text(x(i), y(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end


% Use the plot function to draw the horizontal line
if test == "Miura"
    plot([x(1), x(6)], [y(1), y(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(2), x(5)], [y(2), y(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(5), x(8)], [y(5), y(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(4), x(9)], [y(4), y(9)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    %visualization of the lattice vectors 
    if vis == 1
        plot([x(1) - 0.5*gap, x(1) - 0.5*gap], [y(1), y(7)], 'r-', 'LineWidth', 0.5);
        plot([x(1) - 0.3, x(1) - 0.1], [y(7), y(7)], 'r-', 'LineWidth', 0.5);
        plot([x(1) - 0.3, x(1) - 0.1], [y(1), y(1)], 'r-', 'LineWidth', 0.5);
        h = text(x(1) -0.04, (y(7)+y(1))/1.5, "E2 = " + num2str(y(7) - y(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
        set(h,'Rotation',90);
        
        plot([x(1), x(3)], [y(1)- 0.5*gap, y(1)- 0.5*gap], 'r-', 'LineWidth', 0.5); 
        plot([x(1), x(1)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
        plot([x(3), x(3)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
        text((x(1)+x(3))/1.5, y(1)-0.5*gap, "E1 = " + num2str(x(3) - x(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
    end

end

if test == "Miura4x4" || test == "Miura4x4_diagonal1" || test == "Miura4x4_diagonal2" || test == "Miura4x4_diagonal1_flip" || test == "Miura4x4_diagonal2_flip"
    
    plot([x(1), x(10)], [y(1), y(10)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(2), x(9)], [y(2), y(9)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(3), x(8)], [y(3), y(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(4), x(7)], [y(4), y(7)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(6), x(15)], [y(6), y(15)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(7), x(14)], [y(7), y(14)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(8), x(13)], [y(8), y(13)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(9), x(12)], [y(9), y(12)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(11), x(20)], [y(11), y(20)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(12), x(19)], [y(12), y(19)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(13), x(18)], [y(13), y(18)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(14), x(17)], [y(14), y(17)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(16), x(25)], [y(16), y(25)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(17), x(24)], [y(17), y(24)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(18), x(23)], [y(18), y(23)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(19), x(22)], [y(19), y(22)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    if test == "Miura4x4_diagonal1" || test == "Miura4x4_diagonal2"
        plot([x(1), x(9)], [y(1), y(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(9), x(13)], [y(9), y(13)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(13), x(17)], [y(13), y(17)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(17), x(25)], [y(17), y(25)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

        if test == "Miura4x4_diagonal2"
            plot([x(3), x(7)], [y(3), y(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
            plot([x(7), x(15)], [y(7), y(15)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

            plot([x(11), x(19)], [y(11), y(19)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
            plot([x(19), x(23)], [y(19), y(23)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        end
        
    end

    if test == "Miura4x4_diagonal1_flip"
        plot([x(2), x(10)], [y(2), y(10)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(8), x(12)], [y(8), y(12)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(14), x(18)], [y(14), y(18)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(16), x(24)], [y(16), y(24)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    end

    if test == "Miura4x4_diagonal2_flip"
        plot([x(2), x(10)], [y(2), y(10)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(8), x(12)], [y(8), y(12)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(14), x(18)], [y(14), y(18)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(16), x(24)], [y(16), y(24)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12), x(20)], [y(12), y(20)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(8), x(4)], [y(8), y(4)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(14), x(6)], [y(14), y(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(18), x(22)], [y(18), y(22)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    end

    %visualization of the lattice vectors 
    if vis == 1
        plot([x(1) - 0.5*gap, x(1) - 0.5*gap], [y(1), y(21)], 'r-', 'LineWidth', 0.5);
        plot([x(1) - 0.3, x(1) - 0.1], [y(21), y(21)], 'r-', 'LineWidth', 0.5);
        plot([x(1) - 0.3, x(1) - 0.1], [y(1), y(1)], 'r-', 'LineWidth', 0.5);
        h = text(x(1) -0.04, (y(21)+y(1))/1.5, "E2 = " + num2str(y(21) - y(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
        set(h,'Rotation',90);
        
        plot([x(1), x(5)], [y(1)- 0.5*gap, y(1)- 0.5*gap], 'r-', 'LineWidth', 0.5); 
        plot([x(1), x(1)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
        plot([x(5), x(5)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
        text((x(1)+x(5))/1.5, y(1)-0.5*gap, "E1 = " + num2str(x(5) - x(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
    end
end

if test == "Folding"
    plot([x(1), x(6)], [y(1), y(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(2), x(5)], [y(2), y(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "Planar Kirigami"
    plot([x(2), x(5)], [y(2), y(5)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(5), x(17)], [y(5), y(17)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(2), x(17)], [y(2), y(17)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(6), x(18)], [y(6), y(18)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(6), x(21)], [y(6), y(21)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(7), x(22)], [y(7), y(22)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(7), x(13)], [y(7), y(13)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(23), x(13)], [y(23), y(13)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(3), x(8)], [y(3), y(8)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(3), x(10)], [y(3), y(10)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(17), x(10)], [y(17), y(10)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(17), x(8)], [y(17), y(8)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "Rotating Squares"
    plot([x(1), x(2)], [y(1), y(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(1), x(12)], [y(1), y(12)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(12), x(3)], [y(12), y(3)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(12), x(9)], [y(12), y(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(3), x(6)], [y(3), y(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(7), x(8)], [y(7), y(8)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(6), x(9)], [y(6), y(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    %visualization of the lattice vectors 
    if vis == 1
        plot([x(2) - 0.5*gap, x(2) - 0.5*gap], [y(2), y(10)], 'r-', 'LineWidth', 0.5);
        plot([x(2) - 0.3, x(2) - 0.1], [y(10), y(10)], 'r-', 'LineWidth', 0.5);
        plot([x(2) - 0.3, x(2) - 0.1], [y(2), y(2)], 'r-', 'LineWidth', 0.5);
        h = text(x(2) -0.04, (y(10)+y(2))/1.5, "E2 = " + num2str(y(10) - y(2)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
        set(h,'Rotation',90);
        
        plot([x(1), x(5)], [y(1)- 0.5*gap, y(1)- 0.5*gap], 'r-', 'LineWidth', 0.5); 
        plot([x(1), x(1)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
        plot([x(5), x(5)], [y(1) - 0.3, y(1) - 0.1], 'r-', 'LineWidth', 0.5);
        text((x(1)+x(5))/1.5, y(1)-0.5*gap, "E1 = " + num2str(x(5) - x(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
    end

    % tessellation
    tess = 0; % set to 0 to disable tessellated panels
    if tess
        l1R_x = x(6)-x(12);
        l1R_y = y(6)-y(12);

        plot([x(1)+l1R_x, x(2)+l1R_x], [y(1)+l1R_y, y(2)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l1R_x, x(1)+l1R_x], [y(12)+l1R_y, y(1)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l1R_x, x(3)+l1R_x], [y(12)+l1R_y, y(3)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l1R_x, x(9)+l1R_x], [y(12)+l1R_y, y(9)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(3)+l1R_x, x(6)+l1R_x], [y(3)+l1R_y, y(6)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(6)+l1R_x, x(9)+l1R_x], [y(6)+l1R_y, y(9)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(2)+l1R_x, x(3)+l1R_x], [y(2)+l1R_y, y(3)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(3)+l1R_x, x(4)+l1R_x], [y(3)+l1R_y, y(4)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(4)+l1R_x, x(5)+l1R_x], [y(4)+l1R_y, y(5)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(5)+l1R_x, x(6)+l1R_x], [y(5)+l1R_y, y(6)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(6)+l1R_x, x(7)+l1R_x], [y(6)+l1R_y, y(7)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(7)+l1R_x, x(8)+l1R_x], [y(7)+l1R_y, y(8)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(8)+l1R_x, x(9)+l1R_x], [y(8)+l1R_y, y(9)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(9)+l1R_x, x(10)+l1R_x], [y(9)+l1R_y, y(10)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(10)+l1R_x, x(11)+l1R_x], [y(10)+l1R_y, y(11)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(11)+l1R_x, x(12)+l1R_x], [y(11)+l1R_y, y(12)+l1R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    
        l2R_x = x(11)-x(1);
        l2R_y = y(11)-y(1);
        plot([x(1)+l2R_x, x(2)+l2R_x], [y(1)+l2R_y, y(2)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l2R_x, x(1)+l2R_x], [y(12)+l2R_y, y(1)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l2R_x, x(3)+l2R_x], [y(12)+l2R_y, y(3)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l2R_x, x(9)+l2R_x], [y(12)+l2R_y, y(9)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(3)+l2R_x, x(6)+l2R_x], [y(3)+l2R_y, y(6)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(6)+l2R_x, x(9)+l2R_x], [y(6)+l2R_y, y(9)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(2)+l2R_x, x(3)+l2R_x], [y(2)+l2R_y, y(3)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(3)+l2R_x, x(4)+l2R_x], [y(3)+l2R_y, y(4)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(4)+l2R_x, x(5)+l2R_x], [y(4)+l2R_y, y(5)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(5)+l2R_x, x(6)+l2R_x], [y(5)+l2R_y, y(6)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(6)+l2R_x, x(7)+l2R_x], [y(6)+l2R_y, y(7)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(7)+l2R_x, x(8)+l2R_x], [y(7)+l2R_y, y(8)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(8)+l2R_x, x(9)+l2R_x], [y(8)+l2R_y, y(9)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(9)+l2R_x, x(10)+l2R_x], [y(9)+l2R_y, y(10)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(10)+l2R_x, x(11)+l2R_x], [y(10)+l2R_y, y(11)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(11)+l2R_x, x(12)+l2R_x], [y(11)+l2R_y, y(12)+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

        plot([x(1)+l1R_x+l2R_x, x(2)+l1R_x+l2R_x], [y(1)+l1R_y+l2R_y, y(2)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l1R_x+l2R_x, x(1)+l1R_x+l2R_x], [y(12)+l1R_y+l2R_y, y(1)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l1R_x+l2R_x, x(3)+l1R_x+l2R_x], [y(12)+l1R_y+l2R_y, y(3)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(12)+l1R_x+l2R_x, x(9)+l1R_x+l2R_x], [y(12)+l1R_y+l2R_y, y(9)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(3)+l1R_x+l2R_x, x(6)+l1R_x+l2R_x], [y(3)+l1R_y+l2R_y, y(6)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(6)+l1R_x+l2R_x, x(9)+l1R_x+l2R_x], [y(6)+l1R_y+l2R_y, y(9)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(2)+l1R_x+l2R_x, x(3)+l1R_x+l2R_x], [y(2)+l1R_y+l2R_y, y(3)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(3)+l1R_x+l2R_x, x(4)+l1R_x+l2R_x], [y(3)+l1R_y+l2R_y, y(4)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(4)+l1R_x+l2R_x, x(5)+l1R_x+l2R_x], [y(4)+l1R_y+l2R_y, y(5)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(5)+l1R_x+l2R_x, x(6)+l1R_x+l2R_x], [y(5)+l1R_y+l2R_y, y(6)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(6)+l1R_x+l2R_x, x(7)+l1R_x+l2R_x], [y(6)+l1R_y+l2R_y, y(7)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(7)+l1R_x+l2R_x, x(8)+l1R_x+l2R_x], [y(7)+l1R_y+l2R_y, y(8)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(8)+l1R_x+l2R_x, x(9)+l1R_x+l2R_x], [y(8)+l1R_y+l2R_y, y(9)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(9)+l1R_x+l2R_x, x(10)+l1R_x+l2R_x], [y(9)+l1R_y+l2R_y, y(10)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(10)+l1R_x+l2R_x, x(11)+l1R_x+l2R_x], [y(10)+l1R_y+l2R_y, y(11)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot([x(11)+l1R_x+l2R_x, x(12)+l1R_x+l2R_x], [y(11)+l1R_y+l2R_y, y(12)+l1R_y+l2R_y], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    end
end

if test == "helical_four_triangles"
    plot([x(2), x(5)], [y(2), y(5)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(1), x(3)], [y(1), y(3)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(1), x(4)], [y(1), y(4)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(1), x(5)], [y(1), y(5)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "helical_waterbomb"
    plot([x(1), x(6)], [y(1), y(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(1), x(7)], [y(1), y(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(2), x(7)], [y(2), y(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(3), x(7)], [y(3), y(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(4), x(7)], [y(4), y(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot([x(5), x(7)], [y(5), y(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

for i = 1:length(labels)
    text(x(i), y(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

hold off

end

