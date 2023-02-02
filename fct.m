function out = fct(in)

bool = in < 0;

in(bool) = 0;

out = factorial(in);

out(bool) = NaN;

end