clear;
syms x1 y1 x2 y2 m1 m2

f = 1/(sqrt((m1)^2+(m2)^2));
I = f;
I = int(I, m1, x1, 1+x1);
I = int(I, m2, x2, 1+x2);
I = int(I, y1, 0, 1);
I = int(I, y2, 0, 1);

