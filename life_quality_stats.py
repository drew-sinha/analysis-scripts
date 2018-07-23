# feature = 'centroid_dist'
# cutoff =
# dwell_time = 2  #days

def get_featurespan(worms, feature, cutoff, dwell_time=2, return_crossings=False):

    # TODO: Make the features to negate a little less hardcoded?
    if any([label in feature for label in ['af', 'intensity']]):
        negate = True

    healthspans, crossing_idxs = np.array([]), np.array([])
    for worm in worms:
        ages, feature_data = worm.time_range(
            feature, min_age=0, age_feature='egg_age')
        if negate:
            feature_data *= -1
        adultspan = worm.td.adultspan

        adjusted_data = feature_data - cutoff
        crossings = np.where((adjusted_data[:-1]>0) &
            (adjusted_data[1:] < 0))[0][0]

        if len(crossings) == 0:
            if adjusted_data[0] > 0:
                healthspan = adultspan
                crossing_index = len(adjusted_data)
            elif adjusted_data[0] <= 0:
                healthspan = 0
                crossing_index = 0
        else:
            found_crossing = False
            for crossing_index in crossings:
                dwell_limit = np.where(ages > (ages[crossing_idx] + dwell_time))[0][0]
                if adjusted_data[crossing_index+1:crossing_index+1+dwell_limit] < 0:
                    healthspan = ages[crossing_index]
                    found_crossing = True
                    break
            if not found_crossing:
                healthspan = ages[crossing_index] # Default to the last crossing in a bad situation

        healthspans = np.append(healthspans, healthspan)
        crossing_idxs = np.append(crossing_idxs, crossing_index)

    if return_crossings:
        return healthspans, crossing_idxs
    else:
        return healthspans
