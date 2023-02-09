function diameter = diam_1(len1, len2, distance)
point_1 = [0, len1/2];
point_2 = [0, -len1/2];
point_3 = [distance, len2/2];
%coefficient_a = polyfit(point_1, point_2, 1);
%midpoint_a = 0;
%inter_a = coefficient_a(1);
%slope_a = coefficient_a(2);
%coefficient_b = polyfit(point_1, point_3, 1);
midpoint_b = (point_1 + point_3)/2;
%inter_b = coefficient_b(1);
slope_b = (len2/2 - len1/2)/distance;
%find the perpendicular bisector
%bisector_1_slope = 0;
bisector_2_slope = -1/slope_b;


%eqn2 = -bisector_2_slope*x + y == -bisector_2_slope*midpoint_b(1) + midpoint_b(2);

x = (bisector_2_slope*midpoint_b(1) - midpoint_b(2))/(bisector_2_slope);
y = 0
%find the intercection
diameter = sqrt((x)^2 + (len1/2)^2);

end








