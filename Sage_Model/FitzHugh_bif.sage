reset()

var('x y a b c t')

f1 = diff(x,t) == a - b*x - x^3 + c*y
f2 = diff(y,t) == -y - x

the_field = plot_vector_field((0, f2.subs(a=1,b=1,c=1)), (a, -2, 2), (y, -2, 2))

the_field.show()

