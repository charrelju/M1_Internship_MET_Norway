from scipy.stats.stats import pearsonr

def circular_correlation_coefficient(variable1, variable2):
        #
        var1_rad = variable1 * np.pi / 180
        var2_rad = variable2 * np.pi / 180
        #
        var1_com = np.full(np.shape(variable1), np.nan, dtype = complex)
        var2_com = np.full(np.shape(variable2), np.nan, dtype = complex)
        for i in range(0, len(variable1)):
                var1_com[i] = cmath.rect(1, var1_rad[i])
                var2_com[i] = cmath.rect(1, var2_rad[i])
        #
        var1_com_mean = np.mean(var1_com)
        var2_com_mean = np.mean(var2_com)
        #
        var1_rad_mean = cmath.phase(var1_com_mean)
        if var1_rad_mean < 0:
                var1_rad_mean = var1_rad_mean + 2 * np.pi
        #
        var2_rad_mean = cmath.phase(var2_com_mean)
        if var2_rad_mean < 0:
                var2_rad_mean = var2_rad_mean + 2 * np.pi
        #
        diff_var1 = var1_rad - var1_rad_mean
        diff_var2 = var2_rad - var2_rad_mean
        #
        Rcc = np.sum(np.sin(diff_var1) * np.sin(diff_var2)) / np.sqrt(np.sum(np.sin(diff_var1) ** 2) * np.sum(np.sin(diff_var2) ** 2))
        #
        return Rcc
