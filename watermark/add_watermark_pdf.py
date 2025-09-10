from PyPDF2 import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
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

def add_watermark_to_pdf(input_pdf, output_pdf, watermark_text, owner_password=None):
    """Add watermark to a PDF file"""
    try:
        if not os.path.exists(input_pdf):
            raise FileNotFoundError(f"Input file not found: {input_pdf}")

        reader = PdfReader(input_pdf)
        writer = PdfWriter()

        for page in reader.pages:
            # Create the first watermark (adjusted to fit within the page)
            watermark1 = create_watermark(watermark_text, 200, 400)  # Adjusted y position
            if watermark1:
                page.merge_page(watermark1.pages[0])
                logger.info("First watermark added.")
            else:
                logger.error("Failed to create first watermark.")
            
            # Create the second watermark (adjusted to fit within the page)
            watermark2 = create_watermark(watermark_text, 200, 100)  # Adjusted y position
            if watermark2:
                page.merge_page(watermark2.pages[0])
                logger.info("Second watermark added.")
            else:
                logger.error("Failed to create second watermark.")
                
            writer.add_page(page)

        # Add password protection if owner_password is provided
        if owner_password:
            try:
                # Set permissions to restrict editing but allow reading
                writer.encrypt(
                    user_password="",  # Empty string means no password to open
                    owner_password=owner_password,  # Use the password passed as parameter
                    use_128bit=True,
                    permissions_flag=2048  # Allow printing but restrict editing
                )
                logger.info("PDF encryption added successfully.")
            except Exception as e:
                logger.error(f"Error adding encryption: {e}")
                return False

        with open(output_pdf, "wb") as output_file:
            writer.write(output_file)
        return True
    except Exception as e:
        logger.error(f"Error processing {input_pdf}: {e}")
        return False

def batch_add_watermark(input_folder, output_folder, watermark_text, user_password=None, owner_password=None):
    """Batch process PDF files in a folder and its subfolders"""
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

    for dirpath, _, filenames in os.walk(input_folder):
        for filename in filenames:
            if filename.endswith(".pdf"):
                input_path = os.path.join(dirpath, filename)
                # Create corresponding output path
                relative_path = os.path.relpath(dirpath, input_folder)
                output_path = os.path.join(output_folder, relative_path, filename)

                # Ensure the output subdirectory exists
                os.makedirs(os.path.dirname(output_path), exist_ok=True)

                if add_watermark_to_pdf(input_path, output_path, watermark_text,
                                        owner_password=owner_password):
                    success_count += 1
                    logger.info(f"Successfully processed: {filename}")
                else:
                    failure_count += 1
                    logger.error(f"Failed to process: {filename}")

    logger.info(f"Processing complete. Success: {success_count}, Failures: {failure_count}")
    return failure_count == 0

if __name__ == "__main__":
    input_folder = "/home/yincy/disk14/yr/20250414-山西医科大学-高丽娟/report"
    output_folder = "/home/yincy/disk14/yr/20250414-山西医科大学-高丽娟/report"
    watermark_text = "上海萤锐科技有限公司"
    
    # Only set owner_password to restrict editing, leave user_password as None
    owner_password = "Yr2k5#Bx9m@Pn3v_Kj8wQ$nM4hL7dR2z"  # Password for editing
    user_password = None                    # No password for opening

    if not os.path.exists(input_folder):
        logger.error(f"Input folder does not exist: {input_folder}")
    elif not os.path.isdir(input_folder):
        logger.error(f"Input path is not a directory: {input_folder}")
    else:
        batch_add_watermark(input_folder, output_folder, watermark_text,
                          user_password=user_password,
                          owner_password=owner_password)

