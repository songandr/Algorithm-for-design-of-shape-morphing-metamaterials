function plot1vector3D(vec, name, title1, vis, test)

% plots a column vector of 3-dimensional points (x, y, z1, x2, y2,
% z2....)

% vis is a boolean that determines whether the visualization of the
% symmetry vectors is off or on

gap = 0.4;
vec = cell2mat(vec);
x = zeros(length(vec)/3, 1);
y = x;
z = x;

for i = 1:length(vec)/3
x(i) = vec(3*i-2);
y(i) = vec(3*i-1);
z(i) = vec(3*i);
labels{i} = [name, num2str(i)];
end


figure
xlim([min(x) - gap, max(x) + gap]);
ylim([min(y) - gap, max(y) + gap]);
zlim([min(z) - gap, max(z) + gap]);
plot3(x, y, z);
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');

hold on

for i = 1:length(labels)
    text(x(i), y(i), z(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

title(title1)

%Use the plot function to draw lines

if test == "Miura"
    plot3([x(1), x(6)], [y(1), y(6)], [z(1), z(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(5)], [y(2), y(5)], [z(2), z(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(5), x(8)], [y(5), y(8)], [z(5), z(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(4), x(9)], [y(4), y(9)], [z(4), z(9)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

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
    
    plot3([x(1), x(10)], [y(1), y(10)], [z(1), z(10)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(9)], [y(2), y(9)], [z(2), z(9)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(3), x(8)], [y(3), y(8)], [z(3), z(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(4), x(7)], [y(4), y(7)], [z(4), z(7)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(6), x(15)], [y(6), y(15)], [z(6), z(15)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(7), x(14)], [y(7), y(14)], [z(7), z(14)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(8), x(13)], [y(8), y(13)], [z(8), z(13)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(9), x(12)], [y(9), y(12)], [z(9), z(12)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(11), x(20)], [y(11), y(20)], [z(11), z(20)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12), x(19)], [y(12), y(19)], [z(12), z(19)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(13), x(18)], [y(13), y(18)], [z(13), z(18)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(14), x(17)], [y(14), y(17)], [z(14), z(17)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(16), x(25)], [y(16), y(25)], [z(16), z(25)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(17), x(24)], [y(17), y(24)], [z(17), z(24)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(18), x(23)], [y(18), y(23)], [z(18), z(23)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(19), x(22)], [y(19), y(22)], [z(19), z(22)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    if test == "Miura4x4_diagonal1" || test == "Miura4x4_diagonal2"
        plot3([x(1), x(9)], [y(1), y(9)], [z(1), z(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(9), x(13)], [y(9), y(13)], [z(9), z(13)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(13), x(17)], [y(13), y(17)], [z(13), z(17)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(17), x(25)], [y(17), y(25)], [z(17), z(25)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

        if test == "Miura4x4_diagonal2"
            plot3([x(3), x(7)], [y(3), y(7)], [z(3), z(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
            plot3([x(7), x(15)], [y(7), y(15)], [z(7), z(15)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

            plot3([x(11), x(19)], [y(11), y(19)], [z(11), z(19)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
            plot3([x(19), x(23)], [y(19), y(23)], [z(19), z(23)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        end
        
    end

    if test == "Miura4x4_diagonal1_flip"
        plot3([x(2), x(10)], [y(2), y(10)], [z(2), z(10)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(8), x(12)], [y(8), y(12)], [z(8), z(12)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(14), x(18)], [y(14), y(18)], [z(14), z(18)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(16), x(24)], [y(16), y(24)], [z(16), z(24)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    end

    if test == "Miura4x4_diagonal2_flip"
        plot3([x(2), x(10)], [y(2), y(10)], [z(2), z(10)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(8), x(12)], [y(8), y(12)], [z(8), z(12)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(14), x(18)], [y(14), y(18)], [z(14), z(18)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(16), x(24)], [y(16), y(24)], [z(16), z(24)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(12), x(20)], [y(12), y(20)], [z(12), z(20)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(8), x(4)], [y(8), y(4)], [z(8), z(4)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(14), x(6)], [y(14), y(6)], [z(14), z(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
        plot3([x(18), x(22)], [y(18), y(22)], [z(18), z(22)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
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
    plot3([x(1), x(6)], [y(1), y(6)], [z(1), z(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(5)], [y(2), y(5)], [z(2), z(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "Rotating Squares"
    plot3([x(1), x(2)], [y(1), y(2)], [z(1), z(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12), x(1)], [y(12), y(1)], [z(12), z(1)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12), x(3)], [y(12), y(3)], [z(12), z(3)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12), x(9)], [y(12), y(9)], [z(12), z(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(3), x(6)], [y(3), y(6)], [z(3), z(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(7), x(8)], [y(7), y(8)], [z(7), z(8)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(6), x(9)], [y(6), y(9)], [z(6), z(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

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

    view([0 90]) % top view for a flat y configuration

    % tessellation
    %{
    l1R_x = x(5)-x(1);
    l1R_y = y(5)-y(1);
    plot3([x(1)+l1R_x, x(2)+l1R_x], [y(1)+l1R_y, y(2)+l1R_y], [z(1), z(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12)+l1R_x, x(1)+l1R_x], [y(12)+l1R_y, y(1)+l1R_y], [z(12), z(1)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12)+l1R_x, x(3)+l1R_x], [y(12)+l1R_y, y(3)+l1R_y], [z(12), z(3)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(12)+l1R_x, x(9)+l1R_x], [y(12)+l1R_y, y(9)+l1R_y], [z(12), z(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(3)+l1R_x, x(6)+l1R_x], [y(3)+l1R_y, y(6)+l1R_y], [z(3), z(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(7)+l1R_x, x(8)+l1R_x], [y(7)+l1R_y, y(8)+l1R_y], [z(7), z(8)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(6)+l1R_x, x(9)+l1R_x], [y(6)+l1R_y, y(9)+l1R_y], [z(6), z(9)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    %}
end

if test == "Square Twist"
    plot3([x(1), x(2)], [y(1), y(2)], [z(1), z(2)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(1), x(8)], [y(1), y(8)], [z(1), z(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(7)], [y(2), y(7)], [z(2), z(7)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(7), x(8)], [y(7), y(8)], [z(7), z(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    plot3([x(6), x(7)], [y(6), y(7)], [z(6), z(7)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(3)], [y(2), y(3)], [z(2), z(3)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    plot3([x(3), x(4)], [y(3), y(4)], [z(3), z(4)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(3), x(6)], [y(3), y(6)], [z(3), z(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(4), x(5)], [y(4), y(5)], [z(4), z(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(5), x(6)], [y(5), y(6)], [z(5), z(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    plot3([x(5), x(12)], [y(5), y(12)], [z(5), z(12)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(11), x(6)], [y(11), y(6)], [z(11), z(6)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    
    plot3([x(11), x(12)], [y(11), y(12)], [z(11), z(12)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(11), x(14)], [y(11), y(14)], [z(11), z(14)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(13), x(14)], [y(13), y(14)], [z(13), z(14)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(13), x(12)], [y(13), y(12)], [z(13), z(12)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    plot3([x(11), x(10)], [y(11), y(10)], [z(11), z(10)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(15), x(14)], [y(15), y(14)], [z(15), z(14)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    plot3([x(15), x(10)], [y(15), y(10)], [z(15), z(10)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(15), x(16)], [y(15), y(16)], [z(15), z(16)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);    
    plot3([x(9), x(10)], [y(9), y(10)], [z(9), z(10)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(9), x(16)], [y(9), y(16)], [z(9), z(16)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

    plot3([x(7), x(10)], [y(7), y(10)], [z(7), z(10)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(9), x(8)], [y(9), y(8)], [z(9), z(8)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);

end

if test == "helical_four_triangles"
    plot3([x(2), x(5)], [y(2), y(5)], [z(2), z(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(1), x(3)], [y(1), y(3)], [z(1), z(3)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(1), x(4)], [y(1), y(4)], [z(1), z(4)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(1), x(5)], [y(1), y(5)], [z(1), z(5)],'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

if test == "helical_waterbomb"
    plot3([x(1), x(6)], [y(1), y(6)], [z(1), z(6)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(1), x(7)], [y(1), y(7)], [z(1), z(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(2), x(7)], [y(2), y(7)], [z(2), z(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(3), x(7)], [y(3), y(7)], [z(3), z(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(4), x(7)], [y(4), y(7)], [z(4), z(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    plot3([x(5), x(7)], [y(5), y(7)], [z(5), z(7)], 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
end

%{
for i = 1:length(labels)
    text(x(i), y(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
%}
hold off
end