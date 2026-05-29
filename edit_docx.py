import docx
import re
import sys

try:
    doc_path = r"e:\New folder (8)\UPDATED RESULTS.docx"
    save_path = r"e:\New folder (8)\UPDATED_RESULTS_REVISED.docx"
    doc = docx.Document(doc_path)

    for p in doc.paragraphs:
        # Fix 100-110 J
        if '110 J' in p.text:
            p.text = re.sub(r'100.110\s*J', '100-110 J', p.text)
        
        # Figure 2 change
        orig_fig2 = 'The results demonstrate that material efficiency is a key enabler of the proposed system, as reduced mass directly enhances the ability to convert environmental wind energy into motion. This validates that the structural design choices are aligned with the objective of achieving wind-driven locomotion.'
        new_fig2 = 'Consequently, the achieved material efficiency acts as a primary enabler for converting environmental wind energy into motion, fully validating the structural design choices for wind-driven locomotion.'
        if orig_fig2 in p.text:
            p.text = p.text.replace(orig_fig2, new_fig2)

        # Figure 4 change
        orig_fig4 = 'Such behavior highlights the fundamental operating principle of tensegrity-based motion. This directly validates the design approach, showing that controlled internal actuation can influence global motion without the need for traditional rigid mechanisms.'
        new_fig4 = 'Ultimately, such behavior highlights the fundamental operating principle of tensegrity-based motion, proving that controlled internal actuation can successfully dictate global trajectory without relying on traditional rigid mechanisms.'
        if orig_fig4 in p.text:
            p.text = p.text.replace(orig_fig4, new_fig4)

        # Figure 6 change
        orig_fig6 = 'The ability to reconfigure while maintaining stability highlights a key advantage of tensegrity systems over conventional rigid robots. This demonstrates that the proposed design effectively utilizes structural flexibility as a mechanism for motion generation.'
        new_fig6 = 'This adaptive reconfiguration—achieved without sacrificing overall structural stability—highlights a distinct advantage of tensegrity systems over conventional rigid robots, effectively utilizing structural flexibility as a fundamental mechanism for motion generation.'
        if orig_fig6 in p.text:
            p.text = p.text.replace(orig_fig6, new_fig6)

    doc.save(save_path)
    print(f"Successfully saved revised document to: {save_path}")

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
