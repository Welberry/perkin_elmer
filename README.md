# Perkin Elmer

Codes for processing diffuse X-ray scattering data from a Perkin Elmer detector at the APS Chicago.
Made bad_pixel_map by reading BZC_II_LT/dark_air-02332.dark.tif into ImageJ, 
ran a 5 pixel radius median smoother, reread BZC_II_LT/dark_air-02332.dark.tif
and then did a difference between the smoothed and unsmoothed images.
Arbitrarily decided a difference of 50000+ represented a bad pixel (otherwise
it picks up the sharp lines between the physical sensor blocks that make up the
detector). I subtracted this amount from difference image, changed it to 8-bit
which discards the negative numbers, multiplied by 1000 to make it 0 or 255 and
then inverted it so that the black pixels are the ones to mask. Then saved it
as a pgm file.
