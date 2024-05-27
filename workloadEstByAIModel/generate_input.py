import csv



if __name__ == "__main__":
    # generate csv

    # data set name, number of cycle, dims (3), space at each dim, step size, num of adv steps
    # assuming seeding strategy is global vise
    # should we consider step size?
    # fix step size and number of steps now (2000)
    # do not consider these things in model now
    global_data_set_type=[
        ["astro",1,256,256,256,1.0/256.0,1.0/256.0,1.0/256.0,0.005],
        ["clover",1,256,256,256,1.0/256.0,1.0/256.0,1.0/256.0,0.001],
        ["isabel",1,500,500,100,1,1,1,0.1],
        ["redsea",1,500,500,50,1,1,1,0.1]]
    
    decomp_pattern = [
        [8,2,2,2],
        [16,4,2,2],
        [32,4,4,2],
        [64,4,4,4]
    ]
    
    # title of the csv 
    with open('input.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        attribute = ["dataset_name", "num_sub_blocks_x","num_sub_blocks_y","num_sub_blocks_z","space_x","space_y","space_z","subblock_index_x","subblock_index_y","subblock_index_z","step_size"]
        writer.writerow(attribute)

        # for each data set
        # for each step
        # for each decomposition
        # for id in each decomposed block
        for data_set_info in global_data_set_type:
            data_set_name = data_set_info[0]
            data_set_num = data_set_info[1]
            space_x = data_set_info[5]
            space_y = data_set_info[6]
            space_z = data_set_info[7]
            adv_step_size = data_set_info[8]
            for step in range(data_set_num):
                full_data_set_name=data_set_name+"_"+str(step)
                for decomp_info in decomp_pattern:
                    num_blocks_x=decomp_info[1]
                    num_blocks_y=decomp_info[2]
                    num_blocks_z=decomp_info[3]
                    for index_z in range(num_blocks_z):
                        for index_y in range(num_blocks_y):
                            for index_x in range(num_blocks_x):
                                one_item=[full_data_set_name,num_blocks_x,num_blocks_y,num_blocks_z,space_x,space_y,space_z,index_x,index_y,index_z,adv_step_size]
                                writer.writerow(one_item)
                    


# for each configure, using the actual program to generate a label (which is the workload estimation value)
# issue, different data set has different dims, 
# is it possible to let them generate same number of features based on the 3dconv operations