/* (c) 2014 - 2016 Open Source Geospatial Foundation - all rights reserved
 * (c) 2001 - 2013 OpenPlans
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */
package org.geoserver.wps.gs.download;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.imageio.stream.ImageOutputStream;
import javax.media.jai.BorderExtender;
import javax.media.jai.ImageLayout;
import javax.media.jai.Interpolation;
import javax.media.jai.JAI;
import javax.media.jai.operator.MosaicDescriptor;

import org.geoserver.catalog.CoverageInfo;
import org.geoserver.data.util.CoverageUtils;
import org.geoserver.platform.resource.Resource;
import org.geoserver.wps.WPSException;
import org.geoserver.wps.ppio.ComplexPPIO;
import org.geoserver.wps.ppio.ProcessParameterIO;
import org.geoserver.wps.resource.GridCoverageResource;
import org.geoserver.wps.resource.WPSResourceManager;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridCoverageFactory;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.coverage.grid.io.AbstractGridFormat;
import org.geotools.coverage.grid.io.GridCoverage2DReader;
import org.geotools.coverage.processing.Operations;
import org.geotools.data.Parameter;
import org.geotools.geometry.GeneralEnvelope;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.geotools.process.ProcessException;
import org.geotools.process.raster.CropCoverage;
import org.geotools.referencing.CRS;
import org.geotools.renderer.lite.gridcoverage2d.GridCoverageRenderer;
import org.geotools.resources.coverage.FeatureUtilities;
import org.geotools.styling.RasterSymbolizer;
import org.geotools.styling.RasterSymbolizerImpl;
import org.geotools.util.logging.Logging;
import org.opengis.coverage.grid.GridEnvelope;
import org.opengis.filter.Filter;
import org.opengis.geometry.Envelope;
import org.opengis.parameter.GeneralParameterDescriptor;
import org.opengis.parameter.GeneralParameterValue;
import org.opengis.parameter.ParameterValue;
import org.opengis.parameter.ParameterValueGroup;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.datum.PixelInCell;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;
import org.opengis.util.ProgressListener;
import org.springframework.context.ApplicationContext;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.PrecisionModel;

import it.geosolutions.imageio.stream.output.ImageOutputStreamAdapter;
import it.geosolutions.io.output.adapter.OutputStreamAdapter;

/**
 * Implements the download services for raster data. If limits are configured this class will use {@link LimitedImageOutputStream}, which raises an
 * exception when the download size exceeded the limits.
 * 
 * @author Simone Giannecchini, GeoSolutions SAS
 * 
 */
class RasterDownload {

    private static final Logger LOGGER = Logging.getLogger(RasterDownload.class);

    /** The {@link DownloadServiceConfiguration} object containing the configured limits. */
    private DownloadServiceConfiguration limits;

	/** The resource manager for handling the used resources. */
	private WPSResourceManager resourceManager;

	private final static GridCoverageFactory GC_FACTORY = new GridCoverageFactory();

	private final static RasterSymbolizer RS = new RasterSymbolizerImpl();

    private final static BorderExtender BORDER_EXTENDER_COPY = BorderExtender.createInstance(BorderExtender.BORDER_COPY);

    static final RenderingHints BORDER_EXTENDER_HINTS = new RenderingHints(JAI.KEY_BORDER_EXTENDER, BORDER_EXTENDER_COPY);

    /**
     * The application context used to look-up PPIO factories
     */
    private ApplicationContext context;

    /**
     * Constructor, takes a {@link DownloadEstimatorProcess}.
     * 
     * @param limits the {@link DownloadEstimatorProcess} to check for not exceeding the download
     *        limits.
     * @param resourceManager the {@link WPSResourceManager} to handl generated resources
     * @param context
     */
    public RasterDownload(DownloadServiceConfiguration limits, WPSResourceManager resourceManager,
            ApplicationContext context) {
        this.limits = limits;
        this.resourceManager = resourceManager;
        this.context = context;
    }

    /**
     * This method does the following operations:
     * <ul>
     * <li>Uses only those bands specified by indices (if defined)</li>
     * <li>Reprojection of the coverage (if needed)</li>
     * <li>Clips the coverage (if needed)</li>
     * <li>Scales the coverage to match the target size (if needed)</li>
     * <li>Writes the result</li>
     * <li>Cleanup the generated coverages</li>
     * </ul>
     * 
     * @param mimeType mimetype of the result
     * @param progressListener listener to use for logging the operations
     * @param coverageInfo resource associated to the Coverage
     * @param roi input ROI object
     * @param targetCRS CRS of the file to write
     * @param clip indicates if the clipping geometry must be exactly that of the ROI or simply its envelope
     * @param interpolation interpolation method to use when reprojecting / scaling
     * @param targetSizeX the size of the target image along the X axis
     * @param targetSizeY the size of the target image along the Y axis
     * @param bandIndices the indices of the bands used for the final result
     * @param filter the {@link Filter} to load the data
     *
     */
    public Resource execute(String mimeType, final ProgressListener progressListener,
            CoverageInfo coverageInfo, Geometry roi, CoordinateReferenceSystem targetCRS,
            boolean clip, Filter filter, Interpolation interpolation, Integer targetSizeX,
            Integer targetSizeY, int[] bandIndices) throws Exception {

		GridCoverage2D clippedGridCoverage = null, parentCoverage = null, originalGridCoverage = null,
				reprojectedCoverage = null;
        try {

            //
            // look for output extension. Tiff/tif/geotiff will be all treated as GeoTIFF
            //

            // prepare native CRS
            CoordinateReferenceSystem nativeCRS = DownloadUtilities.getNativeCRS(coverageInfo);
            if (LOGGER.isLoggable(Level.FINE)) {
                LOGGER.fine("Native CRS is " + nativeCRS.toWKT());
            }

            //
            // STEP 0 - Push ROI back to native CRS (if ROI is provided)
            //
            ROIManager roiManager = null;
            if (roi != null) {
                if (LOGGER.isLoggable(Level.FINE)) {
                    LOGGER.log(Level.FINE, "Pushing ROI to native CRS");
                }
                final CoordinateReferenceSystem roiCRS = (CoordinateReferenceSystem) roi
                        .getUserData();
                roiManager = new ROIManager(roi, roiCRS);
            }

            //
            // check reprojection is needed
            //
            boolean reproject = false;
            MathTransform reprojectionTrasform = null;
            if (targetCRS != null && !CRS.equalsIgnoreMetadata(nativeCRS, targetCRS)) {
                if (LOGGER.isLoggable(Level.FINE)) {
                    LOGGER.log(Level.FINE, "Checking if reprojection is needed");
                }
                // testing reprojection...
                reprojectionTrasform = CRS.findMathTransform(nativeCRS, targetCRS, true);
                if (!reprojectionTrasform.isIdentity()) {
                    // avoid doing the transform if this is the identity
                    reproject = true;
                    if (LOGGER.isLoggable(Level.FINE)) {
                        LOGGER.log(Level.FINE, "Reprojection needed");
                    }
                }
            } else {
                targetCRS = nativeCRS;
            }

            // get a reader for this CoverageInfo
            if (LOGGER.isLoggable(Level.FINE)) {
                LOGGER.log(Level.FINE, "Getting reader for the coverage");
            }
            final GridCoverage2DReader reader = (GridCoverage2DReader) coverageInfo
                    .getGridCoverageReader(null, null);
            final ParameterValueGroup readParametersDescriptor = reader.getFormat()
                    .getReadParameters();
            final List<GeneralParameterDescriptor> parameterDescriptors = readParametersDescriptor
                    .getDescriptor().descriptors();
            // get the configured metadata for this coverage without
            Map<String, Serializable> coverageParameters = coverageInfo.getParameters();
            GeneralParameterValue[] readParameters = CoverageUtils.getParameters(
                    readParametersDescriptor, coverageParameters, false);

            // merge support for filter
            if (filter != null) {
                if (LOGGER.isLoggable(Level.FINE)) {
                    LOGGER.log(Level.FINE, "Add the filter");
                }
                readParameters = CoverageUtils.mergeParameter(parameterDescriptors, readParameters,
                        filter, "FILTER", "Filter");
            }

            // read GridGeometry preparation and scaling setup if needed
            GridGeometry2D requestedGridGeometry = null;
            boolean imposedSize = false;
            ReferencedEnvelope targetEnvelope = null;
			if (roi != null) {
				// set crs in roi manager
				roiManager.useNativeCRS(reader.getCoordinateReferenceSystem());
				roiManager.useTargetCRS(targetCRS);
				targetEnvelope = new ReferencedEnvelope(roiManager.getRoiInTargetCRS().getEnvelopeInternal(),
						targetCRS);
				if (LOGGER.isLoggable(Level.FINE)) {
					LOGGER.log(Level.FINE, "Preparing the GridGeometry for cropping input layer with ROI");
				}
				if (targetSizeX == null && targetSizeY == null) {
					// No size is specified. Just do a read and reproject (if needed) + a final crop
					readParameters = computeReadingParams(reader, parameterDescriptors, readParameters, roiManager);
				} else {
					if (targetSizeX == null || targetSizeY == null) {
						// one of the 2 sizes is not specified. Delegate
						// scaleToTarget to compute the second one.
						ScaleToTarget scaling = new ScaleToTarget(reader);
						scaling.setTargetSize(targetSizeX, targetSizeY);
						Integer[] computedSizes = scaling.getTargetSize();
						targetSizeX = computedSizes[0];
						targetSizeY = computedSizes[1];
					}
					imposedSize = true;

					// Since we have imposed a target size, delegate GridCoverageRenderer to do all the job
					requestedGridGeometry = new GridGeometry2D(new GridEnvelope2D(0, 0, targetSizeX, targetSizeY),
							targetEnvelope);
				}
			}

            readParameters = updateReadParams(readParameters, parameterDescriptors, bandIndices);

            // Setting background values and color
            double[] backgroundValues = getBackgroundValues(coverageParameters, readParameters);
            Color color = backgroundValues != null ? new Color((int) backgroundValues[0], (int) backgroundValues[0],
					(int) backgroundValues[0]) : Color.BLACK;

			if (imposedSize) {
				// Delegate the GridCoverageRenderer to do all the work of reprojection/interpolation 
				GridCoverageRenderer renderer = new GridCoverageRenderer(targetCRS, targetEnvelope,
						new Rectangle(0, 0, targetSizeX, targetSizeY), null, null);

				RenderedImage ri = renderer.renderImage(reader, readParameters, RS, interpolation, color, 512, 512);
				if (ri == null) {
					throw new WPSException("The reader did not return anything"
							+ "It normally means there is nothing there, or the data got filtered out by the ROI or filter");
				}

				parentCoverage = (GridCoverage2D) ri.getProperty("ParentCoverage");
			} else {
				// If not, proceed with standard read and reproject
				originalGridCoverage = reader.read(readParameters);

				// check, the reader might have returned a null coverage
				if (originalGridCoverage == null) {
					throw new WPSException("The reader did not return any data for current input "
							+ "parameters. It normally means there is nothing there, or the data got filtered out by the ROI or filter");
				}

				// Reproject if needed
				reprojectedCoverage = originalGridCoverage;
				if (reproject) {
					if (LOGGER.isLoggable(Level.FINE)) {
						LOGGER.log(Level.FINE, "Reprojecting the layer");
					}
					reprojectedCoverage = (GridCoverage2D) Operations.DEFAULT.resample(originalGridCoverage, targetCRS, null,
							interpolation, backgroundValues);
				}
				parentCoverage = reprojectedCoverage;
			}
            // Add a bandSelectProcess call if the reader doesn't support bands
			parentCoverage = checkAndPerformBandSelection(parentCoverage);

            //
            // Clip if needed
            //
            if (roi != null) {
                if (LOGGER.isLoggable(Level.FINE)) {
                    LOGGER.log(Level.FINE, "Cropping the layer");
                }
                // ROI requires a crop/clip
                boolean crop = true;

                if (imposedSize) {
                	crop = false; // There might be the case that GridCoverageRenderer already provided what requested
	                final RenderedImage rasterData = parentCoverage.getRenderedImage();
	                final GridEnvelope requestedRange = (GridEnvelope) requestedGridGeometry.getGridRange();
	
	                // Preliminar check between requested imageLayout and coverage imageLayout
	                final int requestedW = requestedRange.getSpan(0);
	                final int requestedH = requestedRange.getSpan(1);
	                final int imageW = rasterData.getWidth();
	                final int imageH = rasterData.getHeight();

					if ((imageW == requestedW) && (imageH == requestedH) && !clip) {
						// No refining is needed. Write it as is
						return writeRaster(mimeType, coverageInfo, parentCoverage);
						
					} else {
						// Check if an actual crop is needed
						crop = cropIsNeeded(rasterData, requestedRange);
	
						if (!crop && !clip) {
							// The extent of the returned image is smaller than the requested extent. 
							// Let's do a mosaic to return the requested extent instead.  
							clippedGridCoverage = extendToRegion(rasterData, requestedGridGeometry, backgroundValues);
						}
					}
                }

				if (crop || clip) {
					// Crop or Clip
					final CropCoverage cropCoverage = new CropCoverage();

					// Get the proper ROI (depending on clip parameter and CRS)
					Geometry croppingRoi = roiManager.getTargetRoi(clip);
					clippedGridCoverage = cropCoverage.execute(parentCoverage, croppingRoi, progressListener);

					if (clippedGridCoverage == null) {
						throw new WPSException("No data left after applying the ROI. This means there "
								+ "is source data, but none matching the requested ROI");
					}
				}
			} else {
				// do nothing
				clippedGridCoverage = parentCoverage;
			}
            
            //
            // Writing
            //
            return writeRaster(mimeType, coverageInfo, clippedGridCoverage);
        } finally {
            if (originalGridCoverage != null) {
                resourceManager.addResource(new GridCoverageResource(originalGridCoverage));
            }

            if (parentCoverage != null) {
                resourceManager.addResource(new GridCoverageResource(parentCoverage));
            }
            
            if (reprojectedCoverage != null) {
                resourceManager.addResource(new GridCoverageResource(reprojectedCoverage));
            }

            if (clippedGridCoverage != null) {
                resourceManager.addResource(new GridCoverageResource(clippedGridCoverage));
            }

        }
    }

	private double[] getBackgroundValues(Map<String, Serializable> coverageParameters,
			GeneralParameterValue[] readParameters) {
		double[] backgroundValues = null;
		if (coverageParameters != null && coverageParameters.containsKey("BackgroundValues")) {
			for (GeneralParameterValue readParameter : readParameters) {
				if ("BackgroundValues".equalsIgnoreCase(readParameter.getDescriptor().getName().toString())) {
					Object bgValue = ((ParameterValue) readParameter).getValue();
					if (bgValue != null && bgValue instanceof double[]) {
						backgroundValues = ((double[]) bgValue);
					}
					break;
				}
			}
		}
		return backgroundValues;
	}

	private boolean cropIsNeeded(RenderedImage rasterData, GridEnvelope requestedRange) {
		final int requestedW = requestedRange.getSpan(0);
		final int requestedH = requestedRange.getSpan(1);
		final int requestedMinX = requestedRange.getLow(0);
		final int requestedMinY = requestedRange.getLow(1);
		final int requestedMaxX = requestedRange.getHigh(0);
		final int requestedMaxY = requestedRange.getHigh(1);

		final int imageW = rasterData.getWidth();
		final int imageH = rasterData.getHeight();
		final int minX = rasterData.getMinX();
		final int minY = rasterData.getMinY();
		final int maxX = minX + imageW - 1;
		final int maxY = minY + imageH - 1;

		if (((imageW + ((minX >= 0) ? 0 : minX)) >= requestedW) && ((imageH + ((minY >= 0) ? 0 : minY)) >= requestedH)
				&& ((requestedMinX >= minX) && (requestedMinY >= minY) && (requestedMaxX <= maxX)
						&& (requestedMaxY <= maxY))) {
			// The extent of the returned image contains the requested extent.
			// Let's do a crop
			return true;
		}
		return false;
	}

	private GeneralParameterValue[] computeReadingParams(GridCoverage2DReader reader,
			List<GeneralParameterDescriptor> parameterDescriptors, GeneralParameterValue[] readParameters,
			ROIManager roiManager) throws TransformException, IOException {
		final ReferencedEnvelope roiEnvelope = new ReferencedEnvelope(
				roiManager.getSafeRoiInNativeCRS().getEnvelopeInternal(), roiManager.getNativeCRS());
//		final Polygon originalEnvelopeAsPolygon = FeatureUtilities.getPolygon(reader.getOriginalEnvelope(),
//				new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING)));
//		originalEnvelopeAsPolygon.setUserData(roiManager.getNativeCRS());
//		final ReferencedEnvelope originalEnvelope = JTS.toEnvelope(originalEnvelopeAsPolygon);
//		// calculate intersection between original envelope and ROI, as blindly trusting
//		// the ROI may give issues with scaling, if target size is not specified for
//		// both X and Y dimensions
//		final ReferencedEnvelope intersection = originalEnvelope.intersection(roiEnvelope);
//		// take scaling into account
//		ScaleToTarget scaling = new ScaleToTarget(reader, intersection);
		ScaleToTarget scaling = new ScaleToTarget(reader, roiEnvelope);
		scaling.setTargetSize(null, null);
		GridGeometry2D gg2D = scaling.getGridGeometry();
		return CoverageUtils.mergeParameter(parameterDescriptors, readParameters, gg2D,
				AbstractGridFormat.READ_GRIDGEOMETRY2D.getName().getCode());
	}

	private GridCoverage2D extendToRegion(RenderedImage rasterData, GridGeometry2D requestedGridGeometry,
			double[] backgroundValues) {

		final GridEnvelope requestedRange = (GridEnvelope) requestedGridGeometry.getGridRange();
		final int requestedW = requestedRange.getSpan(0);
		final int requestedH = requestedRange.getSpan(1);

		ImageLayout layout = new ImageLayout();
		// Setting output image layout
		// I can impose 0,0 as the origin since we are working under
		// that assumption as part of the reader request contract
		layout.setHeight(requestedH).setWidth(requestedW).setMinX(0).setMinY(0);
		layout.setSampleModel(rasterData.getSampleModel().createCompatibleSampleModel(requestedW, requestedH));
		layout.setColorModel(rasterData.getColorModel());

		final RenderingHints jaiHints = new RenderingHints(JAI.KEY_IMAGE_LAYOUT, layout);
		jaiHints.add(BORDER_EXTENDER_HINTS);
		final RenderedImage mosaic = MosaicDescriptor.create(new RenderedImage[] { rasterData },
				MosaicDescriptor.MOSAIC_TYPE_BLEND, null, null, null, backgroundValues, jaiHints);

		// Setting up a gridCoverage on top of the new image
		final Envelope envelope = new GeneralEnvelope(new GridEnvelope2D(0, 0, requestedW, requestedH),
				PixelInCell.CELL_CENTER, requestedGridGeometry.getGridToCRS(),
				requestedGridGeometry.getCoordinateReferenceSystem());

		return GC_FACTORY.create("mosaiced", mosaic, envelope);

	}

	private GridCoverage2D checkAndPerformBandSelection(GridCoverage2D parentCoverage) {
		// //
		// // STEP 0 - Check for bands, select only those specified
		// //
		// if (bandIndices!=null && bandIndices.length>0){
		// //check band indices are valid
		// int sampleDimensionsNumber =
		// originalGridCoverage.getNumSampleDimensions();
		// for (int i:bandIndices){
		// if (i<0 || i>=sampleDimensionsNumber){
		// throw new WPSException(
		// "Band index "+i+" is invalid for the current input raster. "
		// + "This raster contains "+sampleDimensionsNumber+" band"
		// + (sampleDimensionsNumber>1?"s":""));
		// }
		// }
		// BandSelectProcess bandSelectProcess = new BandSelectProcess();
		//
		// //using null for the VisibleSampleDimension parameter of
		// BandSelectProcess.execute.
		// //GeoTools BandSelector2D takes care of remapping visible band index
		// //or assigns it to first band in order if remapping is not possible
		// bandFilteredCoverage = bandSelectProcess.execute(
		// originalGridCoverage, bandIndices, null);
		//
		// }else{
		// bandFilteredCoverage = originalGridCoverage;
		// }
		return parentCoverage;
	}

	private GeneralParameterValue[] updateReadParams(GeneralParameterValue[] readParameters,
			List<GeneralParameterDescriptor> parameterDescriptors, int[] bandIndices) {

		// make sure we work in streaming fashion
		boolean replacedJai = false;
		for (GeneralParameterValue pv : readParameters) {
			String pdCode = pv.getDescriptor().getName().getCode();
			if (AbstractGridFormat.USE_JAI_IMAGEREAD.getName().getCode().equals(pdCode)) {
				replacedJai = true;
				ParameterValue pvalue = (ParameterValue) pv;
				pvalue.setValue(Boolean.TRUE);
				break;
			}
		}

		if (!replacedJai) {
			readParameters = CoverageUtils.mergeParameter(parameterDescriptors, readParameters, Boolean.TRUE,
					AbstractGridFormat.USE_JAI_IMAGEREAD.getName().getCode());
		}

		// Setting band selection parameter
		boolean replacedBands = false;
		for (GeneralParameterValue pv : readParameters) {
			String pdCode = pv.getDescriptor().getName().getCode();
			if (AbstractGridFormat.BANDS.getName().getCode().equals(pdCode)) {
				replacedBands = true;
				ParameterValue pvalue = (ParameterValue) pv;
				pvalue.setValue(bandIndices);
				break;
			}
		}
		if (!replacedBands) {
			readParameters = CoverageUtils.mergeParameter(parameterDescriptors, readParameters, bandIndices,
					AbstractGridFormat.BANDS.getName().getCode());
		}
		return readParameters;
	}

	/**
     * Writes the providede GridCoverage as a GeoTiff file.
     * 
     * @param mimeType result mimetype
     * @param coverageInfo resource associated to the input coverage
     * @param gridCoverage gridcoverage to write
     * @return a {@link File} that points to the GridCoverage we wrote.
     * 
     */
    private Resource writeRaster(String mimeType, CoverageInfo coverageInfo, GridCoverage2D gridCoverage)
            throws Exception {
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.log(Level.FINE, "Writing raster");
        }
        // limits
        long limit = DownloadServiceConfiguration.NO_LIMIT;
        if (limits.getHardOutputLimit() > 0) {
            limit = limits.getHardOutputLimit();
            if (LOGGER.isLoggable(Level.FINE)) {
                LOGGER.log(Level.FINE, "Hard output limits set to " + limit);
            }
        } else {
            if (LOGGER.isLoggable(Level.FINE)) {
                LOGGER.log(Level.FINE, "Hard output limit unset");
            }
        }

        // Search a proper PPIO
        Parameter<GridCoverage2D> gridParam = new Parameter<GridCoverage2D>("fakeParam",
                GridCoverage2D.class);
        ProcessParameterIO ppio_ = DownloadUtilities.find(gridParam, context, mimeType,
                false);
        if (ppio_ == null) {
            throw new ProcessException("Don't know how to encode in mime type " + mimeType);
        } else if (!(ppio_ instanceof ComplexPPIO)) {
            throw new ProcessException("Invalid PPIO found " + ppio_.getIdentifer());
        }
        final ComplexPPIO complexPPIO = (ComplexPPIO) ppio_;
        String extension = complexPPIO.getFileExtension();

        // writing the output to a temporary folder
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.log(Level.FINE, "Writing file in a temporary folder");
        }
        final Resource output = resourceManager.getTemporaryResource("." + extension);

        // the limit output stream will throw an exception if the process is trying to writer more than the max allowed bytes
        final ImageOutputStream fileImageOutputStreamExtImpl = new ImageOutputStreamAdapter(
                output.out());
        ImageOutputStream os = null;
        // write
        try {
            // If limit is defined, LimitedImageOutputStream is used
            if (limit > DownloadServiceConfiguration.NO_LIMIT) {
                os = new LimitedImageOutputStream(fileImageOutputStreamExtImpl, limit) {

                    @Override
                    protected void raiseError(long pSizeMax, long pCount) throws IOException {
                        IOException e = new IOException(
                                "Download Exceeded the maximum HARD allowed size!");
                        throw e;
                    }
                };
            } else {
                os = fileImageOutputStreamExtImpl;
            }
            // Encoding the GridCoverage
            complexPPIO.encode(gridCoverage, new OutputStreamAdapter(os));
            os.flush();
        } finally {
            try {
                if (os != null) {
                    os.close();
                }
            } catch (Exception e) {
                if (LOGGER.isLoggable(Level.FINE)) {
                    LOGGER.log(Level.FINE, e.getLocalizedMessage(), e);
                }
            }
        }
        return output;
    }
}