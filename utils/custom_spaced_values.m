% Generate a skewed grid with denser sampling near smaller values.
function lambda_values = custom_spaced_values(startvalue, endvalue, num, exponent)
    % Generate custom spaced values skewed towards smaller values
    values = linspace(startvalue^(1/exponent), endvalue^(1/exponent), num) .^ exponent;
    lambda_values = values;
end
