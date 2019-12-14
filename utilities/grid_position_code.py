def acquire_grid_positions(scope, spacing_cc, grid_dims, cycle_through_positions=True):
    '''
        spacing_cc: spacing center to center between positions in the same row
        grid_dims: [num_x_positions, num_y_positions]
        cycle_through_positions: bool flag to specify whether to cycle through positions and set z manually;
            if False, defaults to z of the top-left position specified by user
    '''

    # TODO: Figure out how to invert positions so that one can generate the shortest path across the grid

    try:
        input('Specify top-left position in the grid')
        tl_position = scope.stage.position

        input('Specify the bottom right position in the grid')
        br_position = scope.stage.position

        grid_x = np.linspace(0, spacing_cc*grid_dims[0],grid_dims[0]+1)
        grid_y = np.linspace(0, spacing_cc*grid_dims[1],grid_dims[1]+1)

        xcoor, ycoor = np.meshgrid(grid_x, grid_y)
        grid_points = numpy.array([position for position in zip(xcoor.flatten(), ycoor.flatten())])
        tl_grid_position, br_grid_position = grid_points[0], grid_points[-1]

        _, _, _, grid_points = fit_affine.fit_affine([tl_grid_position, br_grid_position], [tl_position, br_position], find_translation=True, find_scale=False) # Need to map grid points to user-provided positions; should allow reflection?
        assert grid_points.min() > 0

        if cycle_through_positions:
            positions = []
            print('Cycling through calculated positions.')
            for num, (x, y) in enumerate(grid_points):
                scope.stage.position = [x, y, scope.stage.z]
                input('')
                positions.append(scope.stage.position)
                print(f'Position {num}: {scope.stage.positions[-1]}', end='')
        else:
            positions = [[x, y, scope.stage.position] for (x,y) in grid_points]

    except KeyboardInterrupt:
        break

    return positions
