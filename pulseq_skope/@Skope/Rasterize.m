% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function t = Rasterize(this, t, type)
% shift value 't' to raster-grid defined by 'type'
    
    if (type == "grad") % gradient raster time
        raster_time = this.sys.gradRasterTime;
    else
        error("Raster-type unknown!");
    end
    
    t =  (t / raster_time) * raster_time;
end 