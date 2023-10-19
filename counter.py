if __name__ == '__main__':
    
    number_of_struct    = 0
    number_of_iteration = 0

    with open('mp-20_out.cif','r') as mp:

        for line in mp:

            number_of_iteration += 1
            
            if line.startswith('data_'):
                number_of_struct += 1

            if number_of_iteration % 1000 == 0:
               print(F"Number of lines read so far: {number_of_iteration}")

        print(f"Total number of lines in mp-20_out.cif: {number_of_iteration}")
        print(f"Total number of structures in mp-20_out.cif : {number_of_struct}")
