YOUNGS_MODULUS = 70.0e9
POISSON_RATIO = 0.33
DENSITY = 2700.0
YIELD_STRESS = 250.0e6
ULTIMATE_STRESS = 300.0e6


def shear_modulus(youngs_modulus=YOUNGS_MODULUS, poisson_ratio=POISSON_RATIO):
    return youngs_modulus / (2.0 * (1.0 + poisson_ratio))


def flexural_rigidity(thickness, youngs_modulus=YOUNGS_MODULUS, poisson_ratio=POISSON_RATIO):
    return youngs_modulus * thickness**3 / (12.0 * (1.0 - poisson_ratio**2))
