# Thermal Radiation Diffusion Input Parser



module Inputs

function readInputs(filename)
        """ Function to read the input file and store the parameters in a dictionary
            Parameters:
            filename: String - Path to the input file
            Returns:
            params: Dict - Dictionary of input parameters
        """
        print("Reading input file... \n")

        #filename = raw"src\inputs\test_input.txt"

        params = Dict()

        current_key = ""
        current_value = []
        in_array = false

        open(filename, "r") do f
            for line in eachline(f)
                line = strip(line)  # Remove any surrounding whitespace
                
                if isempty(line) || startswith(line, "#")  # Skip empty lines and comments
                    continue
                end
    
                # Check for a new key-value pair (e.g., "key = value")
                if !in_array && occursin(r"^\w+ = ", line)
                    # Split the key and value
                    parts = split(line, "=")
                    current_key = strip(parts[1])
                    current_value = [strip(parts[2])]
                    in_array = false
                    # Save the current value if we have one
                    if !isempty(current_key) && !startswith(current_value[1], "[")
                        params[current_key] = current_value[1]
                    elseif startswith(current_value[1], "[")  # Start of an array
                        if !endswith(strip(current_value[1]), "]")  # Check if array ends this line
                            in_array = true
                        end
                        old_value = current_value[1]
                        current_value = []
                        push!(current_value, strip(old_value, ['[', ']']))  # Remove the opening "[" and closing "]" if it exists before store the value
                        if in_array == false
                            params[current_key] = current_value # Save if array is on one line
                        end
                    end
                elseif in_array  # Array elements spread across lines
                    if endswith(strip(line), "]")  # End of an array
                        push!(current_value, strip(line, ']'))  # Remove the closing "]"
                        params[current_key] = current_value  # Save the key-value pair
                        in_array = false
                    else
                        push!(current_value, strip(line))
                    end
                end
            end
        end

         # Check for separated array quantities and join them into a single array
        for (key, value) in params
                if contains(join(value), ',')
                    if key == raw"MESHNODES" || key == raw"XMESHNODES" || key == raw"YMESHNODES"
                        params[key] = separated_array_parser(value) # Needed to accurately calculate dx/dy values for fine non-uniform meshes
                    else
                        params[key] = separated_array_parser(value)
                    end
                end
        end

    return params
end

   function separated_array_parser(input_array)
        """ Function to parse an array of values separated by commas
            Parameters:
            input_array: Array{String} - Array of values separated by commas
            Returns:
            output_array: Array{precision} - Array of parsed values
        """
        if contains(join(input_array), '(')
            output_array = [eval(Meta.parse(join(input_array)))]
        else
            output_array = parse.(Float64, split(join(input_array), ","))
        end
        return output_array
    end



end