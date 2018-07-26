import numpy as np

def get_featurespan(worms, feature, cutoff, dwell_time=2, return_crossings=False):
    '''Get amount of time feature is remains above/below a prescribed cutoff value;
        assumes that all ages are in *days*.

        Parameters
            worms - Worms object of worm data
            feature - str denoting feature of interest in the worms object
            cutoff - float value for cutoff (in units of the feature)
            dwell_time - (optional) minimum time an individual must stay below the
                cutoff to continue in "featurespan" (in days)
            return_crossings - optional bool denoting whether to return the feature value
                at the time which an individual passed into poor health (good for
                debugging and visualization)

        Returns
            numpy array of featurespans (in days)
            (optional, if return_crossings is True) crossing values

    '''

    # Ensure that feature decreases with time
    # TODO: Make the features to negate a little less hardcoded?
    if any([label in feature for label in ['af', 'intensity']]):
        negate = True
    else: negate = False

    featurespans, crossing_vals = np.array([]), np.array([])
    for worm in worms:
        ages, feature_data = worm.get_time_range(feature, min_age=0, age_feature='adult_age')
        if negate:
            feature_data *= -1

        adjusted_data = feature_data - cutoff
        crossings = np.where((adjusted_data[:-1] > 0) & (adjusted_data[1:] < 0))[0]

        if len(crossings) == 0:
            if adjusted_data[0] > 0:
                featurespan = worm.adultspan
                crossing_val = feature_data[-1]
            elif adjusted_data[0] <= 0:
                featurespan = 0
                crossing_val = feature_data[0]
        else:
            found_crossing = False
            for crossing_index in crossings:
                dwell_end = np.where(ages > (ages[crossing_index] + dwell_time))[0][0]
                if (adjusted_data[crossing_index+1:dwell_end] < 0).all():
                    featurespan = ages[crossing_index]
                    crossing_val = feature_data[crossing_index]
                    found_crossing = True
                    break
            if not found_crossing:
                featurespan = ages[crossing_index] # Default to the last crossing in a bad situation
                crossing_val = feature_data[crossing_index]

        featurespans = np.append(featurespans, featurespan)
        crossing_vals = np.append(crossing_vals, crossing_val)

    if return_crossings:
        return featurespans, crossing_vals
    else:
        return featurespans
