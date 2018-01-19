reset()

var('x y a b c t')
c = 1
a = 1
f1 = diff(x,t) == a - x - x^3 + c*y
f2 = diff(y,t) == -y - x

plot_vector_field((0, f2), (a, -2, 2), (y, -2, 2))

