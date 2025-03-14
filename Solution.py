import pandas as pd
import json
from math import atan, sin, cos, tan, sqrt, radians, atan2

# We can use Haversine formula or Vincenty formula to calculate the distance between two points on the Earth's surface.
# While Haversine formula assume the Earth is a perfect sphere [https://en.wikipedia.org/wiki/Haversine_formula]
# Vincenty formula are based on the assumption that Earth is an oblate spheroid [https://en.wikipedia.org/wiki/Vincenty%27s_formulae]
# I have chosen Vincenty formula due to its higher accuracy compared to Haversine formula.

def vincenty_distance(lat1, lon1, lat2, lon2):
    # WGS-84 constants provided by Wikipedia
    a = 6378137.0  # Semi-major axis (meters)
    f = 1 / 298.257223563  # Flattening
    b = (1 - f) * a  # Length of semi-minor axis

    # Convert lat and lon into radians
    lat1, lon1, lat2, lon2 = radians(lat1), radians(lon1), radians(lat2), radians(lon2)

    # U1 and U2 should be the reduced latitudes
    U1 = atan((1 - f) * tan(lat1))
    U2 = atan((1 - f) * tan(lat2))

    L = lon2 - lon1
    Lambda = L 

    # Iteration settings
    MAX_ITERATIONS = 1000
    CONVERGENCE_THRESHOLD = 1e-12  # Stop when change is small

    # Calculate U1, U2 and L and set the initial value of Lambda = L then evaluate the following equations until Lambda converges
    for _ in range(MAX_ITERATIONS):
        sin_lambda = sin(Lambda)
        cos_lambda = cos(Lambda)

        sin_sigma = sqrt((cos(U2) * sin_lambda) ** 2 + (cos(U1) * sin(U2) - sin(U1)* cos(U2)*cos_lambda) ** 2)
        cos_sigma = sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos_lambda
        sigma = atan2(sin_sigma, cos_sigma)
        sin_alpha = cos(U1) * cos(U2) * sin_lambda / sin_sigma
        cossquared_alpha = 1 - sin_alpha ** 2
        cos_2sigma_m = cos_sigma - 2 * sin(U1) * sin(U2) / cossquared_alpha
        C = f / 16 * cossquared_alpha * (4 + f * (4 - 3 * cossquared_alpha))
        Lambda_prev = Lambda
        Lambda = L + (1 - C) * f * sin_alpha * (sigma + C * sin_sigma * (cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m ** 2)))
        
        # If lambda converges to desired degree of accuracy of 1e-12, break the loop
        if abs(Lambda - Lambda_prev) < CONVERGENCE_THRESHOLD:
            break
    else:
        return None
        
    u_squared = cossquared_alpha * (a ** 2 - b ** 2) / b ** 2
    A = 1 + u_squared / 16384 * (4096 + u_squared * (-768 + u_squared * (320 - 175 * u_squared)))
    B = u_squared / 1024 * (256 + u_squared * (-128 + u_squared * (74 - 47 * u_squared)))
    delta_sigma = B * sin_sigma * (cos_2sigma_m + B / 4 * (cos_sigma * (-1 + 2 * cos_2sigma_m ** 2) - B / 6 * 
                                                           cos_2sigma_m * (-3 + 4 * sin_sigma ** 2) * (-3 + 4 * 
                                                                                                       cos_2sigma_m ** 2)))
    # This is the distance 
    s = b * A * (sigma - delta_sigma)
    return s

# That took a while...

def solution():
    sensor1_data = pd.read_csv('SensorData1.csv')

    with open('SensorData2.json') as f:
        sensor2_data = json.load(f)

    sensor2_df = pd.DataFrame(sensor2_data)

    matched_pairs = {}

    for _, row1 in sensor1_data.iterrows():
        for _, row2 in sensor2_df.iterrows():
            distance = vincenty_distance(row1['latitude'], row1['longitude'], row2['latitude'], row2['longitude'])
            if distance <= 100:
                matched_pairs[row1['id']] = row2['id']

    with open('SolutionOutput.json', 'w') as f:
        json.dump(matched_pairs, f, indent= 4)
    
solution()