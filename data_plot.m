data = load('data_position.mat');
data = data.data;



experiment_name = data { 4 }. name ;
y = data { 4 }. y ;
u = data { 4 }. u ;




figure;
plot(u);
title('Signal u');

figure;
plot(y);
title('Signal y');