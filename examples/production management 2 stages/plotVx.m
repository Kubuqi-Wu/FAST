% Plot value function + objective at first node
function plotVx

[x1, x2] = meshgrid(0:0.1:5,0:0.1:5) ;
d = [2 3 4 5] ;
p = 1/4 ;
z = zeros(size(x1)) ;
obj = zeros(size(x1)) ;
for i = 1:size(x1, 1)
    for j = 1:size(x1, 2)
        obj(i,j) = x1(i,j) + 2 * x2(i,j) ;
        z(i,j) = 0 ;
        for k = 1:length(d)
            s2 = min(x2(i,j), d(k)) ;
            s1 = min(x1(i,j), d(k) - s2) ;
            z(i,j) = z(i,j) + p * (- 2 * s1 - 3 * s2) ;
        end
    end
end

figure ; 
surf(x1, x2, z + obj) ;
xlabel('x1') ;
ylabel('x2') ;
zlabel('z') ;