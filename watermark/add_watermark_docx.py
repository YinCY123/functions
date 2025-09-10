from PyPDF2 import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from docx import Document
from docx.shared import Pt, RGBColor
from docx.oxml.ns import qn
from docx.oxml import parse_xml
import io
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Chinese font configuration
CHINESE_FONT_PATHS = [
    "/usr/share/fonts/truetype/arphic/uming.ttc",  # Linux
    "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",  # macOS
    "C:/Windows/Fonts/simsun.ttc",  # Windows
    "simsun.ttc"  # Fallback
]

def register_chinese_font():
    """Register a Chinese font for PDF generation"""
    for path in CHINESE_FONT_PATHS:
        if os.path.exists(path):
            try:
                pdfmetrics.registerFont(TTFont('ChineseFont', path))
                return True
            except Exception as e:
                logger.warning(f"Failed to register font {path}: {e}")
    logger.error("No suitable Chinese font found")
    return False

def create_watermark(watermark_text, x, y, rotation_angle=45):
    """Create a watermark PDF with Chinese text support"""
    try:
        packet = io.BytesIO()
        can = canvas.Canvas(packet, pagesize=letter)
        
        # Use Chinese font if available, fallback to Helvetica
        if register_chinese_font():
            can.setFont('ChineseFont', 30)
        else:
            can.setFont("Helvetica", 30)
            logger.warning("Using fallback font, Chinese characters may not display correctly")
            
        can.setFillColorRGB(1, 0, 0, alpha=0.3)
        
        # Save state, rotate, draw text, then restore state
        can.saveState()
        can.translate(x, y)
        can.rotate(rotation_angle)
        can.drawString(0, 0, watermark_text)
        can.restoreState()
        can.save()
        packet.seek(0)
        return PdfReader(packet)
    except Exception as e:
        logger.error(f"Error creating watermark PDF: {e}")
        return None

def add_watermark_to_pdf(input_pdf, output_pdf, watermark_text):
    """Add watermark to a PDF file"""
    try:
        if not os.path.exists(input_pdf):
            raise FileNotFoundError(f"Input file not found: {input_pdf}")

        reader = PdfReader(input_pdf)
        writer = PdfWriter()

        for page in reader.pages:
            watermark = create_watermark(watermark_text, 100, 400)
            if watermark:
                page.merge_page(watermark.pages[0])
            writer.add_page(page)

        with open(output_pdf, "wb") as output_file:
            writer.write(output_file)
        return True
    except Exception as e:
        logger.error(f"Error processing {input_pdf}: {e}")
        return False

def batch_add_watermark(input_folder, output_folder, watermark_text):
    """Batch process PDF files in a folder"""
    if not os.path.exists(input_folder):
        logger.error(f"Input folder does not exist: {input_folder}")
        return False

    if not os.path.exists(output_folder):
        try:
            os.makedirs(output_folder)
        except Exception as e:
            logger.error(f"Error creating output folder: {e}")
            return False

    success_count = 0
    failure_count = 0

    for filename in os.listdir(input_folder):
        if filename.endswith(".pdf"):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)
            
            if add_watermark_to_pdf(input_path, output_path, watermark_text):
                success_count += 1
                logger.info(f"Successfully processed: {filename}")
            else:
                failure_count += 1
                logger.error(f"Failed to process: {filename}")

    logger.info(f"Processing complete. Success: {success_count}, Failures: {failure_count}")
    return failure_count == 0

if __name__ == "__main__":
    input_folder = "/home/yincy/disk14/qs/20250210"  # Replace with actual path
    output_folder = "/home/yincy/disk14/qs/20250210"  # Replace with actual path
    watermark_text = "上海萤锐科技有限公司"  # Chinese watermark text

    if not os.path.exists(input_folder):
        logger.error(f"Input folder does not exist: {input_folder}")
    elif not os.path.isdir(input_folder):
        logger.error(f"Input path is not a directory: {input_folder}")
    else:
        batch_add_watermark(input_folder, output_folder, watermark_text)

