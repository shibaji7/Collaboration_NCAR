import numpy

class CalcFov(object):
    """
    Class to calculate fov coords!
    This is mostly copied from DaViTPy.
    """
    
    def __init__(self, frang=180.0, rsep=45.0, hdw=None, nbeams=None,
                 ngates=None, bmsep=None, recrise=None, siteLat=None,
                 siteLon=None, siteBore=None, siteAlt=None, siteYear=None,
                 elevation=None, altitude=300., hop=None, model='IS',
                 date_time=None, coord_alt=0., fov_dir='front'):
        # Define class constants
        rn = 'fov'

        # Test that we have enough input arguments to work with
        if(not hdw and None in [nbeams, ngates, bmsep, recrise, siteLat,
                                 siteLon, siteBore, siteAlt, siteYear]):
            estr = '{:s}: must provide either a hdw object or '.format(rn)
            estr = '{:s}[nbeams, ngates, bmsep, recrise, siteLat,'.format(estr)
            estr = '{:s} siteLon, siteBore, siteAlt, siteYear].'.format(estr)
            logging.error(estr)
            return

        # date_time checking is handled by coord_conv, and it already
        # knows all of the possible coord systems, so no need to do it
        # here.

        # Then assign variables from the hdw object if necessary
        if hdw:
            if not nbeams:
                nbeams = hdw.beams
            if not ngates:
                ngates = hdw.gates
            if not bmsep:
                bmsep = hdw.beam_seperation
            if not recrise:
                recrise = hdw.rx_rise_time
            if not siteLat:
                siteLat = hdw.geographic.lat
            if not siteLon:
                siteLon = hdw.geographic.lon
            if not siteAlt:
                siteAlt = hdw.geographic.alt
            if not siteBore:
                siteBore = hdw.boresight
#             if not siteYear:
#                 siteYear = hdw.time.yr

        # Some type checking is neccessary. If frang, rsep or recrise are
        # arrays, then they should be of shape (nbeams,).  Set a flag if any of
        # frang, rsep or recrise is an array
        is_param_array = False
        if isinstance(frang, numpy.ndarray):
            is_param_array = True
            # Array is adjusted to add on extra beam edge by copying the last
            # element
            if len(frang) != nbeams:
                estr = "{:s}: frang must be a scalar or numpy ".format(rn)
                estr = "{:s}ndarray of size (nbeams). Using first".format(estr)
                estr = "{:s} element: {}".format(estr, frang[0])
                logging.error(estr)
                frang = frang[0] * numpy.ones(nbeams + 1)
            else:
                frang = numpy.append(frang, frang[-1])
        else:
            frang = numpy.array([frang])
        if isinstance(rsep, numpy.ndarray):
            is_param_array = True
            # Array is adjusted to add on extra beam edge by copying the last
            # element
            if len(rsep) != nbeams:
                estr = "{:s}: rsep must be a scalar or numpy ndarray".format(
                    rn)
                estr = "{:s} of size (nbeams). Using first element".format(
                    estr)
                estr = "{:s}: {}".format(estr, rsep[0])
                logging.error(estr)
                rsep = rsep[0] * numpy.ones(nbeams + 1)
            else:
                rsep = numpy.append(rsep, rsep[-1])
        else:
            rsep = numpy.array([rsep])
        if isinstance(recrise, numpy.ndarray):
            is_param_array = True
            # Array is adjusted to add on extra beam edge by copying the last
            # element
            if len(recrise) != nbeams:
                estr = "{:s}: recrise must be a scalar or numpy ".format(rn)
                estr = "{:s}ndarray of size (nbeams). Using first ".format(
                    estr)
                estr = "{:s}element: {}".format(estr, recrise[0])
                logging.error(estr)
                recrise = recrise[0] * numpy.ones(nbeams + 1)
            else:
                recrise = numpy.append(recrise, recrise[-1])
        else:
            recrise = numpy.array([recrise])

        # If altitude, elevation, or hop are arrays, then they should be of
        # shape (nbeams, ngates)
        if isinstance(altitude, numpy.ndarray):
            if altitude.ndim == 1:
                # Array is adjusted to add on extra beam/gate edge by copying
                # the last element and replicating the whole array as many
                # times as beams
                if altitude.size != ngates:
                    estr = '{:s}: altitude must be of a scalar or '.format(rn)
                    estr = '{:s}numpy ndarray of size (ngates) or '.format(
                        estr)
                    estr = '{:s}(nbeans,ngates). Using first '.format(estr)
                    estr = '{:s}element: {}'.format(estr, altitude[0])
                    logging.error(estr)
                    altitude = altitude[0] * numpy.ones((nbeams + 1, ngates + 1))
                else:
                    altitude = numpy.resize(numpy.append(altitude, altitude[-1]),
                                         (nbeams + 1, ngates + 1))
            elif altitude.ndim == 2:
                # Array is adjusted to add on extra beam/gate edge by copying
                # the last row and column
                if altitude.shape != (nbeams, ngates):
                    estr = '{:s}: altitude must be of a scalar or '.format(rn)
                    estr = '{:s}numpy ndarray of size (ngates) or '.format(
                        estr)
                    estr = '{:s}(nbeans,ngates). Using first '.format(estr)
                    estr = '{:s}element: {}'.format(altitude[0])
                    logging.error(estr)
                    altitude = altitude[0] * numpy.ones((nbeams + 1, ngates + 1))
                else:
                    altitude = numpy.append(altitude,
                                         altitude[-1, :].reshape(1, ngates),
                                         axis=0)
                    altitude = numpy.append(altitude,
                                         altitude[:, -1].reshape(nbeams, 1),
                                         axis=1)
            else:
                estr = '{:s}: altitude must be of a scalar or '.format(rn)
                estr = '{:s}numpy ndarray of size (ngates) or '.format(estr)
                estr = '{:s}(nbeans,ngates). Using first element: '.format(
                    estr)
                estr = '{:s}{}'.format(estr, altitude[0])
                logging.error(estr)
                altitude = altitude[0] * numpy.ones((nbeams + 1, ngates + 1))
        if isinstance(elevation, numpy.ndarray):
            if elevation.ndim == 1:
                # Array is adjusted to add on extra beam/gate edge by copying
                # the last element and replicating the whole array as many
                # times as beams
                if elevation.size != ngates:
                    estr = '{:s}: elevation must be of a scalar or '.format(rn)
                    estr = '{:s}numpy ndarray of size (ngates) or '.format(
                        estr)
                    estr = '{:s}(nbeans,ngates). Using first '.format(estr)
                    estr = '{:s}element: {}'.format(estr, elevation[0])
                    logging.error(estr)
                    elevation = elevation[0] * \
                        numpy.ones((nbeams + 1, ngates + 1))
                else:
                    elevation = numpy.resize(numpy.append(elevation, elevation[-1]),
                                          (nbeams + 1, ngates + 1))
            elif elevation.ndim == 2:
                # Array is adjusted to add on extra beam/gate edge by copying
                # the last row and column
                if elevation.shape != (nbeams, ngates):
                    estr = '{:s}: elevation must be of a scalar or '.format(rn)
                    estr = '{:s}numpy ndarray of size (ngates) or '.format(
                        estr)
                    estr = '{:s}(nbeans,ngates). Using first '.format(estr)
                    estr = '{:s}element: {}'.format(estr, elevation[0])
                    logging.error(estr)
                    elevation = elevation[0] * \
                        numpy.ones((nbeams + 1, ngates + 1))
                else:
                    elevation = numpy.append(elevation,
                                          elevation[-1, :].reshape(1, ngates),
                                          axis=0)
                    elevation = numpy.append(elevation,
                                          elevation[:, -1].reshape(nbeams, 1),
                                          axis=1)
            else:
                estr = '{:s}: elevation must be a scalar or '.format(rn)
                estr = '{:s}numpy ndarray of size (ngates) or '.format(estr)
                estr = '{:s}(nbeans,ngates). Using first element'.format(estr)
                estr = '{:s}: {}'.format(estr, elevation[0])
                logging.error(estr)
                elevation = elevation[0] * numpy.ones((nbeams + 1, ngates + 1))

        if isinstance(hop, numpy.ndarray):
            if hop.ndim == 1:
                # Array is adjusted to add on extra beam/gate edge by copying
                # the last element and replicating the whole array as many
                # times as beams
                if hop.size != ngates:
                    estr = '{:s}: hop must be of a scalar or numpy '.format(rn)
                    estr = '{:s}ndarray of size (ngates) or '.format(estr)
                    estr = '{:s}(nbeans,ngates). Using first '.format(estr)
                    estr = '{:s}element: {}'.format(estr, hop[0])
                    logging.error(estr)
                    hop = hop[0] * numpy.ones((nbeams + 1, ngates + 1))
                else:
                    hop = numpy.resize(numpy.append(hop, hop[-1]),
                                    (nbeams + 1, ngates + 1))
            elif hop.ndim == 2:
                # Array is adjusted to add on extra beam/gate edge by copying
                # the last row and column
                if hop.shape != (nbeams, ngates):
                    estr = '{:s}: hop must be of a scalar or numpy '.format(rn)
                    estr = '{:s}ndarray of size (ngates) or '.format(estr)
                    estr = '{:s}(nbeans,ngates). Using first '.format(estr)
                    estr = '{:s}element: {}'.format(hop[0])
                    logging.error(estr)
                    hop = hop[0] * numpy.ones((nbeams + 1, ngates + 1))
                else:
                    hop = numpy.append(hop, hop[-1, :].reshape(1, ngates), axis=0)
                    hop = numpy.append(hop, hop[:, -1].reshape(nbeams, 1), axis=1)
            else:
                estr = '{:s}: hop must be a scalar or numpy ndarray'.format(rn)
                estr = '{:s} of size (ngates) or (nbeams,ngates).'.format(estr)
                estr = '{:s} Using first element: {}'.format(estr, hop[0])
                logging.error(estr)
                hop = hop[0] * numpy.ones((nbeams + 1, ngates + 1))

        # Do for coord_alt what we just did for altitude.
        if isinstance(coord_alt, numpy.ndarray):
            if coord_alt.ndim == 1:
                # Array is adjusted to add on extra beam/gate edge by copying
                # the last element and replicating the whole array as many
                # times as beams
                if coord_alt.size != ngates:
                    estr = '{:s}: coord_alt must be a scalar or '.format(rn)
                    estr = '{:s}numpy ndarray of size (ngates) or '.format(
                        estr)
                    estr = '{:s}(nbeans,ngates). Using first '.format(estr)
                    estr = '{:s}element: {}'.format(estr, coord_alt[0])
                    logging.error(estr)
                    coord_alt = coord_alt[0] * \
                        numpy.ones((nbeams + 1, ngates + 1))
                else:
                    coord_alt = numpy.resize(numpy.append(coord_alt, coord_alt[-1]),
                                          (nbeams + 1, ngates + 1))
            elif coord_alt.ndim == 2:
                # Array is adjusted to add on extra beam/gate edge by copying
                # the last row and column
                if coord_alt.shape != (nbeams, ngates):
                    estr = '{:s}: coord_alt must be a scalar or '.format(estr)
                    estr = '{:s}numpy ndarray of size (ngates) or '.format(
                        estr)
                    estr = '{:s}(nbeans,ngates). Using first '.format(estr)
                    estr = '{:s}element: {}'.format(estr, coord_alt[0])
                    logging.error(estr)
                    coord_alt = coord_alt[0] * \
                        numpy.ones((nbeams + 1, ngates + 1))
                else:
                    coord_alt = numpy.append(coord_alt,
                                          coord_alt[-1, :].reshape(1, ngates),
                                          axis=0)
                    coord_alt = numpy.append(coord_alt,
                                          coord_alt[:, -1].reshape(nbeams, 1),
                                          axis=1)
            else:
                estr = '{:s}: coord_alt must be a scalar or '.format(rn)
                estr = '{:s}numpy ndarray of size (ngates) or '.format(estr)
                estr = '{:s}(nbeans,ngates). Using first element'.format(estr)
                estr = '{:s}: {}'.format(estr, coord_alt[0])
                logging.error(estr)
                coord_alt = coord_alt[0] * numpy.ones((nbeams + 1, ngates + 1))

        # Generate beam/gate arrays
        beams = numpy.arange(nbeams + 1)
        gates = numpy.arange(ngates + 1)

        # Create output arrays
        slant_range_full = numpy.zeros((nbeams + 1, ngates + 1), dtype='float')
        lat_full = numpy.zeros((nbeams + 1, ngates + 1), dtype='float')
        lon_full = numpy.zeros((nbeams + 1, ngates + 1), dtype='float')
        slant_range_center = numpy.zeros((nbeams + 1, ngates + 1), dtype='float')
        lat_center = numpy.zeros((nbeams + 1, ngates + 1), dtype='float')
        lon_center = numpy.zeros((nbeams + 1, ngates + 1), dtype='float')

        # Calculate deviation from boresight for center of beam
        boff_center = bmsep * (beams - (nbeams - 1) / 2.0)
        # Calculate deviation from boresight for edge of beam
        boff_edge = bmsep * (beams - (nbeams - 1) / 2.0 - 0.5)

        # Iterates through beams
        for ib in beams:
            # if none of frang, rsep or recrise are arrays, then only execute
            # this for the first loop, otherwise, repeat for every beam
            if (~is_param_array and ib == 0) or is_param_array:
                # Calculate center slant range
                srang_center = self.slantRange(frang[ib], rsep[ib], recrise[ib],
                                          gates, center=True)
                # Calculate edges slant range
                srang_edge = self.slantRange(frang[ib], rsep[ib], recrise[ib],
                                        gates, center=False)
            # Save into output arrays
            slant_range_center[ib, :-1] = srang_center[:-1]
            slant_range_full[ib, :] = srang_edge

            # Calculate coordinates for Edge and Center of the current beam
            for ig in gates:
                # Handle array-or-not question.
                talt = altitude[ib, ig] if isinstance(altitude, numpy.ndarray) \
                    else altitude
                telv = elevation[ib, ig] if isinstance(elevation, numpy.ndarray) \
                    else elevation
                t_c_alt = coord_alt[ib, ig] \
                    if isinstance(coord_alt, numpy.ndarray) else coord_alt
                thop = hop[ib, ig] if isinstance(hop, numpy.ndarray) else hop

                if model == 'GS':
                    if (~is_param_array and ib == 0) or is_param_array:
                        slant_range_center[ib, ig] = \
                            self.gsMapSlantRange(srang_center[ig], altitude=None,
                                            elevation=None)
                        slant_range_full[ib, ig] = \
                            self.gsMapSlantRange(srang_edge[ig], altitude=None,
                                            elevation=None)
                        srang_center[ig] = slant_range_center[ib, ig]
                        srang_edge[ig] = slant_range_full[ib, ig]

                if (srang_center[ig] != -1) and (srang_edge[ig] != -1):
                    # Then calculate projections
                    latc, lonc = self.calcFieldPnt(siteLat, siteLon, siteAlt * 1e-3,
                                              siteBore, boff_center[ib],
                                              srang_center[ig], elevation=telv,
                                              altitude=talt, hop=thop,
                                              model=model, fov_dir=fov_dir)
                    late, lone = self.calcFieldPnt(siteLat, siteLon, siteAlt * 1e-3,
                                              siteBore, boff_edge[ib],
                                              srang_edge[ig], elevation=telv,
                                              altitude=talt, hop=thop,
                                              model=model, fov_dir=fov_dir)
                else:
                    latc, lonc = numpy.nan, numpy.nan
                    late, lone = numpy.nan, numpy.nan

                # Save into output arrays
                lat_center[ib, ig] = latc
                lon_center[ib, ig] = lonc
                lat_full[ib, ig] = late
                lon_full[ib, ig] = lone

        # Output is...
        self.latCenter = lat_center[:-1, :-1]
        self.lonCenter = lon_center[:-1, :-1]
        self.slantRCenter = slant_range_center[:-1, :-1]
        self.latFull = lat_full
        self.lonFull = lon_full
        self.slantRFull = slant_range_full
        self.beams = beams[:-1]
        self.gates = gates[:-1]
        self.coords = "geo"
        self.fov_dir = fov_dir
        self.model = model

    # *************************************************************
    def __str__(self):
        outstring = 'latCenter: {}\nlonCenter: {}\nlatFull: {}\nlonFull: {} \
                     \nslantRCenter: {}\nslantRFull: {}\nbeams: {} \
                     \ngates: {} \ncoords: {} \nfield of view: {}\
                     \nmodel: {}'.format(numpy.shape(self.latCenter),
                                         numpy.shape(self.lonCenter),
                                         numpy.shape(self.latFull),
                                         numpy.shape(self.lonFull),
                                         numpy.shape(self.slantRCenter),
                                         numpy.shape(self.slantRFull),
                                         numpy.shape(self.beams),
                                         numpy.shape(self.gates), self.coords,
                                         self.fov_dir, self.model)
        return outstring            

    def gsMapSlantRange(self, slant_range, altitude=None, elevation=None):
        """Calculate the ground scatter mapped slant range.
        See Bristow et al. [1994] for more details. (Needs full reference)
        Parameters
        ----------
        slant_range
            normal slant range [km]
        altitude : Optional[float]
            altitude [km] (defaults to 300 km)
        elevation : Optional[float]
            elevation angle [degree]
        Returns
        -------
        gsSlantRange
            ground scatter mapped slant range [km] (typically slightly less than
            0.5 * slant_range.  Will return -1 if
            (slant_range**2 / 4. - altitude**2) >= 0. This occurs when the scatter
            is too close and this model breaks down.

        Shameless Ripoff from DaViTPy
        """
        Re = 6371.0

        # Make sure you have altitude, because these 2 projection models rely on it
        if not elevation and not altitude:
            # Set default altitude to 300 km
            altitude = 300.0
        elif elevation and not altitude:
            # If you have elevation but not altitude, then you calculate altitude,
            # and elevation will be adjusted anyway
            altitude = numpy.sqrt(Re ** 2 + slant_range ** 2 + 2. * slant_range * Re *
                               numpy.sin(numpy.radians(elevation))) - Re

        if (slant_range**2) / 4. - altitude ** 2 >= 0:
            gsSlantRange = Re * \
                numpy.arcsin(numpy.sqrt(slant_range ** 2 / 4. - altitude ** 2) / Re)
            # From Bristow et al. [1994]
        else:
            gsSlantRange = -1

        return gsSlantRange


    def slantRange(self, frang, rsep, recrise, range_gate, center=True):
        """ Calculate slant range
        Parameters
        ----------
        frang : (float)
            first range gate position [km]
        rsep : (float)
            range gate separation [km]
        recrise : (float)
            receiver rise time [us]
        range_gate : (int)
            range gate number(s)
        center : (bool)
            whether or not to compute the slant range in the center of
            the gate rather than at the edge (default=True)
        Returns
        -------
        srang : (float)
            slant range [km]

        Shameless Ripoff from DaViTPy
        """
        # Lag to first range gate [us]
        lagfr = frang * 2.0 / 0.3
        # Sample separation [us]
        smsep = rsep * 2.0 / 0.3
        # Range offset if calculating slant range at center of the gate
        range_offset = -0.5 * rsep if not center else 0.0

        # Slant range [km]
        srang = (lagfr - recrise + range_gate * smsep) * 0.3 / 2.0 + range_offset

        return srang


    def calcFieldPnt(self, tr_glat, tr_glon, tr_alt, boresight, beam_off, slant_range,
                     adjusted_sr=True, elevation=None, altitude=None, hop=None,
                     model=None, coords='geo', gs_loc="G", max_vh=400.0,
                     fov_dir='front', eval_loc=False):
        """Calculate coordinates of field point given the radar coordinates and
        boresight, the pointing direction deviation from boresight and elevation
        angle, and the field point slant range and altitude. Either the elevation
        or the altitude must be provided. If none is provided, the altitude is set
        to 300 km and the elevation evaluated to accomodate altitude and range.
        Parameters
        ----------
        tr_glat
            transmitter latitude [degree, N]
        tr_glon
            transmitter longitude [degree, E]
        tr_alt
            transmitter altitude [km]
        boresight
            boresight azimuth [degree, E]
        beam_off
            beam azimuthal offset from boresight [degree]
        slant_range
            slant range [km]
        adjusted_sr : Optional(bool)
            Denotes whether or not the slant range is the total measured slant
            range (False) or if it has been adjusted to be the slant distance to
            last ionospheric reflection point (True).  (default=True)
        elevation : Optional[float]
            elevation angle [degree] (estimated if None)
        altitude : Optional[float]
            altitude [km] (default 300 km)
        hop : Optional[float]
            backscatter hop (ie 0.5, 1.5 for ionospheric; 1.0, 2.0 for ground)
        model : Optional[str]
            IS : for standard ionopsheric scatter projection model (ignores hop)
            GS : for standard ground scatter projection model (ignores hop)
            S : for standard projection model (uses hop)
            E1 : for Chisham E-region 1/2-hop ionospheric projection model
            F1 : for Chisham F-region 1/2-hop ionospheric projection model
            F3 : for Chisham F-region 1-1/2-hop ionospheric projection model
            C : for Chisham projection model (ionospheric only, ignores hop,
                requires total measured slant range)
            None : if you trust your elevation or altitude values. more to come
        coords
            'geo' (more to come)
        gs_loc : (str)
            Provide last ground scatter location 'G' or ionospheric refraction
            location 'I' for groundscatter (default='G')
        max_vh : (float)
            Maximum height for longer slant ranges in Standard model (default=400)
        fov_dir : (str)
            'front' (default) or 'back'.  Specifies fov direction
        eval_loc : (bool)
            Evaluate the calcualted location based on reasonable tolerances (True)
            or accept the first calculation (False).  Using True gives better
            locations, but restricts data at the furthest range gates.
            (default=False)
        Returns
        ---------
        geo_dict['geoLat'] : (float or numpy.ndarray)
            Field point latitude(s) in degrees or numpy.nan if error
        geo_dict['geoLon'] : (float or numpy.ndarray)
            Field point longitude(s) in degrees or numpy.nan if error

        Shameless Ripoff from DaViTPy
        """
    #     import sys
    #     sys.path.append('../utils/')
        import geoPack
        import model_vheight as vhm

        # Only geo is implemented.
        if coords != "geo":
            logging.error("Only geographic (geo) is implemented in calcFieldPnt.")
            return numpy.nan, numpy.nan

        # Use model to get altitude if desired
        xalt = numpy.nan
        calt = None
        if model is not None:
            if model in ['IS', 'GS', 'S']:
                # The standard model can be used with or without an input altitude
                # or elevation.  Returns an altitude that has been adjusted to
                # comply with common scatter distributions
                if hop is None:
                    if model == "S":
                        # Default to ionospheric backscatter if hop not specified
                        hop = 0.5
                    else:
                        hop = 0.5 if model == "IS" else 1.0

                xalt = vhm.standard_vhm(slant_range, adjusted_sr=adjusted_sr,
                                        max_vh=max_vh, hop=hop, alt=altitude,
                                        elv=elevation)
            else:
                # The Chisham model uses only the total slant range to determine
                # altitude based on years of backscatter data at SAS.  Performs
                # significantly better than the standard model for ionospheric
                # backscatter, not tested on groundscatter
                if adjusted_sr:
                    logging.error("Chisham model needs total slant range")
                    return numpy.nan, numpy.nan

                # Use Chisham model to calculate virtual height
                cmodel = None if model == "C" else model
                xalt, shop = vhm.chisham_vhm(slant_range, cmodel, hop_output=True)

                # If hop is not known, set using model divisions
                if hop is None:
                    hop = shop

                # If hop is greater than 1/2, the elevation angle needs to be
                # calculated from the ground range rather than the virtual height
                if hop > 0.5:
                    calt = float(xalt)

        elif elevation is None or numpy.isnan(elevation):
            if hop is None or adjusted_sr:
                logging.error("Total slant range and hop needed with measurements")
                return numpy.nan, numpy.nan
            if altitude is None or numpy.isnan(altitude):
                logging.error("No observations supplied")
                return numpy.nan, numpy.nan

            # Adjust slant range if there is groundscatter and the location
            # desired is the ionospheric reflection point
            asr = slant_range
            if hop == numpy.floor(hop) and gs_loc == "I":
                asr *= 1.0 - 1.0 / (2.0 * hop)

            # Adjust altitude if it's unrealistic
            if asr < altitude:
                altitude = asr - 10

            xalt = altitude

        # Use model altitude to determine elevation angle and then the location,
        # or if elevation angle was supplied, find the location
        if not numpy.isnan(xalt):
            # Since we have a modeled or measured altitude, start by setting the
            # Earth radius below field point to Earth radius at radar
            (lat, lon, tr_rad) = geoPack.geodToGeoc(tr_glat, tr_glon)
            rad_pos = tr_rad

            # Iterate until the altitude corresponding to the calculated elevation
            # matches the desired altitude.  Assumes straight-line path to last
            # ionospheric scattering point, so adjust slant range if necessary
            # for groundscatter
            asr = slant_range
            shop = hop
            if not adjusted_sr and hop == numpy.floor(hop) and gs_loc == "I":
                asr *= 1.0 - 1.0 / (2.0 * hop)
                shop = hop - 0.5

            # Set safty counter and iteratively determine location
            maxn = 30
            hdel = 100.0
            htol = 0.5
            if (slant_range >= 800.0 and model != 'GS') or shop > 1.0:
                htol = 5.0
            n = 0
            while n < maxn:
                tr_dist = tr_rad + tr_alt
                if calt is not None:
                    # Adjust elevation angle for any hop > 1 (Chisham et al. 2008)
                    pos_dist = rad_pos + calt
                    phi = numpy.arccos((tr_dist**2 + pos_dist**2 - asr**2) /
                                    (2.0 * tr_dist * pos_dist))
                    beta = numpy.arcsin((tr_dist * numpy.sin(phi / (shop * 2.0))) /
                                     (asr / (shop * 2.0)))
                    tel = numpy.pi / 2.0 - beta - phi / (shop * 2.0)

                    if xalt == calt:
                        xalt = numpy.sqrt(tr_rad**2 + asr**2 + 2.0 * asr * tr_rad *
                                       numpy.sin(tel)) - tr_rad
                    tel = numpy.degrees(tel)
                else:
                    # pointing elevation (spherical Earth value) [degree]
                    tel = numpy.arcsin(((rad_pos + xalt)**2 - tr_dist**2 - asr**2) /
                                    (2.0 * tr_dist * asr))
                    tel = numpy.degrees(tel)

                # estimate off-array-normal azimuth (because it varies slightly
                # with elevation) [degree]
                boff = self.calcAzOffBore(tel, beam_off, fov_dir=fov_dir)
                # pointing azimuth
                taz = boresight + boff
                # calculate position of field point
                geo_dict = geoPack.calcDistPnt(tr_glat, tr_glon, tr_alt,
                                               dist=asr, el=tel, az=taz)

                # Update Earth radius
                rad_pos = geo_dict['distRe']

                # stop if the altitude is what we want it to be (or close enough)
                new_hdel = abs(xalt - geo_dict['distAlt'])
                if new_hdel <= htol or not eval_loc:
                    break

                # stop unsuccessfully if the altitude difference hasn't improved
                if abs(new_hdel - hdel) < 1.0e-3:
                    n = maxn

                # Prepare the next iteration
                hdel = new_hdel
                n += 1

            if n >= maxn:
                estr = 'Accuracy on height calculation ({}) not '.format(htol)
                estr = '{:s}reached quick enough. Returning nan, nan.'.format(estr)
                logging.warning(estr)
                return numpy.nan, numpy.nan
            else:
                return geo_dict['distLat'], geo_dict['distLon']
        elif elevation is not None:
            # No projection model (i.e., the elevation or altitude is so good that
            # it gives you the proper projection by simple geometric
            # considerations). Using no models simply means tracing based on
            # trustworthy elevation or altitude
            if hop is None or adjusted_sr:
                logging.error("Hop and total slant range needed with measurements")
                return numpy.nan, numpy.nan

            if numpy.isnan(elevation):
                logging.error("No observations provided")
                return numpy.nan, numpy.nan

            shop = hop - 0.5 if hop == numpy.floor(hop) and gs_loc == "I" else hop
            asr = slant_range
            if hop > 0.5 and hop != shop:
                asr *= 1.0 - 1.0 / (2.0 * hop)

            # The tracing is done by calcDistPnt
            boff = self.calcAzOffBore(elevation, beam_off, fov_dir=fov_dir)
            geo_dict = geoPack.calcDistPnt(tr_glat, tr_glon, tr_alt, dist=asr,
                                           el=elevation, az=boresight + boff)

            return geo_dict['distLat'], geo_dict['distLon']


    def calcAzOffBore(self, elevation, boff_zero, fov_dir='front'):
        """Calculate off-boresight azimuth as a function of elevation angle and
        zero-elevation off-boresight azimuth.
        See Milan et al. [1997] for more details on how this works.
        Parameters
        ----------
        elevation
            elevation angle [degree]
        boff_zero
            zero-elevation off-boresight azimuth [degree]
        fov_dir
            field-of-view direction ('front','back'). Default='front'
        Returns
        -------
        bore_offset
            off-boresight azimuth [degree]
        """
        # Test to see where the true beam direction lies
        bdir = numpy.cos(numpy.radians(boff_zero))**2 - numpy.sin(numpy.radians(elevation))**2

        # Calculate the front fov azimuthal angle off the boresite
        if bdir < 0.0:
            bore_offset = numpy.pi / 2.
        else:
            tan_boff = numpy.sqrt(numpy.sin(numpy.radians(boff_zero))**2 / bdir)
            bore_offset = numpy.arctan(tan_boff)

    # Old version
    #   if bdir < 0.0:
    #       if boff_zero >= 0: boreOffset = numpy.pi/2.
    #       else: boreOffset = -numpy.pi/2.
    #   else:
    #       tan_boff = numpy.sqrt(numpy.sin(numpy.radians(boff_zero))**2 / bdir)
    #       if boff_zero >= 0: boreOffset = atan(tan_boff)
    #       else: boreOffset = -atan(tan_boff)

        # If the rear lobe is desired, adjust the azimuthal offset from the
        # boresite
        if fov_dir is 'back':
            bore_offset = numpy.pi - bore_offset

        # Correct the sign based on the sign of the zero-elevation off-boresight
        # azimuth
        if boff_zero < 0.0:
            bore_offset *= -1.0

        return numpy.degrees(bore_offset)
