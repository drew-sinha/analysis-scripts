import numpy as np

from zplib import fit_affine

def acquire_grid_positions(scope, spacing_cc, grid_dims, cycle_through_positions=True):
    '''
        spacing_cc: spacing center to center between positions in the same row
        grid_dims: [num_x_positions, num_y_positions]
        cycle_through_positions: bool flag to specify whether to cycle through positions and set z manually;
            if False, defaults to z of the top-left position specified by user
            
        Good spacings 
            7x divots: 2.713 (should redo this on a fresh corral)
    '''
    # TODO: Debug issues with wrong direction when using Zach's code.
    # TODO: Figure out how to invert positions so that one can generate the shortest path across the grid

    try:
        input('Specify top-left position in the grid')
        tl_position = scope.stage.position[:-1]

        input('Specify the bottom right position in the grid')
        br_position = scope.stage.position[:-1]

        grid_x = np.linspace(0, spacing_cc*(grid_dims[0]-1),grid_dims[0])
        grid_y = np.linspace(0, spacing_cc*(grid_dims[1]-1),grid_dims[1])

        xcoor, ycoor = np.meshgrid(grid_x, grid_y)
        grid_points = np.array([position for position in zip(xcoor.flatten(), ycoor.flatten())])
        tl_grid_position, br_grid_position = grid_points[0], grid_points[-1]
        #grid_points = np.array([points if ((point_num % (grid_dims[0] * 2)) < grid_dims[0]) else points[::-1] for point_num, points in enumerate(grid_points)]) # Switch every other row to get the shortest path through the grid.

        trans = np.array(tl_position)
        disp = np.array(br_position) - np.array(tl_position)
        theta = np.math.atan2(br_grid_position[1], br_grid_position[0]) - np.math.atan2(disp[1], disp[0])
        rot = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]]) # *think* the signage is right
        grid_points = np.dot(grid_points, rot)  + trans
        
        # Rotations with Zach's code; for now leave it out.
        #rot, _, trans, _ = fit_affine.fit_affine([tl_grid_position, br_grid_position], [tl_position, br_position], find_translation=True, find_scale=False, allow_reflection=True) # Need to map grid points to user-provided positions
        #new_grid_points = np.dot(grid_points, rot) + trans


        if cycle_through_positions:
            positions = []
            print('Cycling through calculated positions. Press enter when position aligned. Press "s"+enter to skip a position.')
            num_acquired_pos = 0
            for num, (x, y) in enumerate(grid_points):
                scope.stage.position = [x, y, scope.stage.z]
                resp = input('')
                if resp == 's':
                    print('Skipping position', end='')
                else:
                    positions.append(scope.stage.position)
                    num_acquired_pos += 1
                    print(f'Position {num_acquired_pos}: {positions[-1]}', end='')
        else:
            positions = [[x, y, scope.stage.position] for (x,y) in grid_points]

    except KeyboardInterrupt:
        return

    return positions
