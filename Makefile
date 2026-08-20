# ===================================
# STRS Project Makefile
# Figure Generation Only
# ===================================

# R executable
R := Rscript

# Directories
FIGURE_DIR := figures
SCRIPT_DIR := scripts/figures

# Figure targets
FIGURES := $(FIGURE_DIR)/Figure1.pdf \
           $(FIGURE_DIR)/Figure2.pdf \
           $(FIGURE_DIR)/Figure3.pdf \
           $(FIGURE_DIR)/Figure4.pdf \
           $(FIGURE_DIR)/Figure5.pdf \
           $(FIGURE_DIR)/FigureS1.pdf \
           $(FIGURE_DIR)/FigureS2.pdf \
           $(FIGURE_DIR)/FigureS4.pdf \
           $(FIGURE_DIR)/FigureS5.pdf \
           $(FIGURE_DIR)/FigureS6.pdf \
           $(FIGURE_DIR)/FigureS8.pdf \
           $(FIGURE_DIR)/FigureS9.pdf \
           $(FIGURE_DIR)/FigureS10.pdf \
           $(FIGURE_DIR)/FigureS11.pdf \
           $(FIGURE_DIR)/FigureS12.pdf \
           $(FIGURE_DIR)/FigureS13.pdf

# ===================================
# MAIN TARGETS
# ===================================

# Default: Generate all figures
.PHONY: all
all: $(FIGURES)
	@echo "All figures generated successfully!"

# Main figures
.PHONY: main
main: $(FIGURE_DIR)/Figure1.pdf \
      $(FIGURE_DIR)/Figure2.pdf \
      $(FIGURE_DIR)/Figure3.pdf \
      $(FIGURE_DIR)/Figure4.pdf \
      $(FIGURE_DIR)/Figure5.pdf
	@echo "Main figures (1-5) generated successfully!"

# Supplementary figures
.PHONY: supplementary
supplementary: $(FIGURE_DIR)/FigureS1.pdf \
               $(FIGURE_DIR)/FigureS2.pdf \
               $(FIGURE_DIR)/FigureS4.pdf \
               $(FIGURE_DIR)/FigureS5.pdf \
               $(FIGURE_DIR)/FigureS6.pdf \
               $(FIGURE_DIR)/FigureS8.pdf \
               $(FIGURE_DIR)/FigureS9.pdf \
               $(FIGURE_DIR)/FigureS10.pdf \
               $(FIGURE_DIR)/FigureS11.pdf \
               $(FIGURE_DIR)/FigureS12.pdf \
               $(FIGURE_DIR)/FigureS13.pdf
	@echo "Supplementary figures generated successfully!"

# ===================================
# FIGURE GENERATION RULES
# ===================================

$(FIGURE_DIR)/Figure1.pdf: $(SCRIPT_DIR)/Figure1.R
	@echo "Generating Figure 1..."
	@$(R) $<

$(FIGURE_DIR)/Figure2.pdf: $(SCRIPT_DIR)/Figure2.R
	@echo "Generating Figure 2..."
	@$(R) $<

$(FIGURE_DIR)/Figure3.pdf: $(SCRIPT_DIR)/Figure3.R
	@echo "Generating Figure 3..."
	@$(R) $<

$(FIGURE_DIR)/Figure4.pdf: $(SCRIPT_DIR)/Figure4.R
	@echo "Generating Figure 4..."
	@$(R) $<

$(FIGURE_DIR)/Figure5.pdf: $(SCRIPT_DIR)/Figure5.R
	@echo "Generating Figure 5..."
	@$(R) $<

$(FIGURE_DIR)/FigureS1.pdf: $(SCRIPT_DIR)/FigureS1.R
	@echo "Generating Supplementary Figure 1..."
	@$(R) $<

$(FIGURE_DIR)/FigureS2.pdf: $(SCRIPT_DIR)/FigureS2.R
	@echo "Generating Supplementary Figure 2..."
	@$(R) $<

$(FIGURE_DIR)/FigureS4.pdf: $(SCRIPT_DIR)/FigureS4.R
	@echo "Generating Supplementary Figure 4..."
	@$(R) $<

$(FIGURE_DIR)/FigureS5.pdf: $(SCRIPT_DIR)/FigureS5.R
	@echo "Generating Supplementary Figure 5..."
	@$(R) $<

$(FIGURE_DIR)/FigureS6.pdf: $(SCRIPT_DIR)/FigureS6.R
	@echo "Generating Supplementary Figure 6..."
	@$(R) $<

$(FIGURE_DIR)/FigureS8.pdf: $(SCRIPT_DIR)/FigureS8.R
	@echo "Generating Supplementary Figure 8..."
	@$(R) $<

$(FIGURE_DIR)/FigureS9.pdf: $(SCRIPT_DIR)/FigureS9.R
	@echo "Generating Supplementary Figure 9..."
	@$(R) $<

$(FIGURE_DIR)/FigureS10.pdf: $(SCRIPT_DIR)/FigureS10.R
	@echo "Generating Supplementary Figure 10..."
	@$(R) $<

$(FIGURE_DIR)/FigureS11.pdf: $(SCRIPT_DIR)/FigureS11.R
	@echo "Generating Supplementary Figure 11..."
	@$(R) $<

$(FIGURE_DIR)/FigureS12.pdf: $(SCRIPT_DIR)/FigureS12.R
	@echo "Generating Supplementary Figure 12..."
	@$(R) $<

$(FIGURE_DIR)/FigureS13.pdf: $(SCRIPT_DIR)/FigureS13.R
	@echo "Generating Supplementary Figure 13..."
	@$(R) $<

# ===================================
# INDIVIDUAL FIGURE SHORTCUTS
# ===================================

.PHONY: fig1 fig2 fig3 fig4 fig5 figs1 figs2 figs4 figs5 figs6 figs8 figs9 figs10 figs11 figs12 figs13
fig1: $(FIGURE_DIR)/Figure1.pdf
fig2: $(FIGURE_DIR)/Figure2.pdf
fig3: $(FIGURE_DIR)/Figure3.pdf
fig4: $(FIGURE_DIR)/Figure4.pdf
fig5: $(FIGURE_DIR)/Figure5.pdf
figs1: $(FIGURE_DIR)/FigureS1.pdf
figs2: $(FIGURE_DIR)/FigureS2.pdf
figs4: $(FIGURE_DIR)/FigureS4.pdf
figs5: $(FIGURE_DIR)/FigureS5.pdf
figs6: $(FIGURE_DIR)/FigureS6.pdf
figs8: $(FIGURE_DIR)/FigureS8.pdf
figs9: $(FIGURE_DIR)/FigureS9.pdf
figs10: $(FIGURE_DIR)/FigureS10.pdf
figs11: $(FIGURE_DIR)/FigureS11.pdf
figs12: $(FIGURE_DIR)/FigureS12.pdf
figs13: $(FIGURE_DIR)/FigureS13.pdf

# ===================================
# MAINTENANCE
# ===================================

.PHONY: clean
clean:
	@echo "Removing all generated figures..."
	@rm -f $(FIGURES)
	@echo "Done!"

.PHONY: clean-plots
clean-plots:
	@echo "Removing intermediate plots directory..."
	@rm -rf plots/
	@echo "Done!"

# ===================================
# HELP
# ===================================

.PHONY: help
help:
	@echo "STRS Project Makefile - Figure Generation"
	@echo "=========================================="
	@echo ""
	@echo "Available commands:"
	@echo "  make              - Generate all figures (default)"
	@echo "  make main         - Generate main figures (Figure 1-5)"
	@echo "  make supplementary - Generate supplementary figures (S1, S2, S4-S6, S8-S13)"
	@echo "  make fig1         - Generate Figure 1"
	@echo "  make fig2         - Generate Figure 2"
	@echo "  make fig3         - Generate Figure 3"
	@echo "  make fig4         - Generate Figure 4"
	@echo "  make fig5         - Generate Figure 5"
	@echo "  make figs1        - Generate Supplementary Figure 1"
	@echo "  make figs2        - Generate Supplementary Figure 2"
	@echo "  make figs4        - Generate Supplementary Figure 4"
	@echo "  make figs5        - Generate Supplementary Figure 5"
	@echo "  make figs6        - Generate Supplementary Figure 6"
	@echo "  make figs8        - Generate Supplementary Figure 8"
	@echo "  make figs9        - Generate Supplementary Figure 9"
	@echo "  make figs10       - Generate Supplementary Figure 10"
	@echo "  make figs11       - Generate Supplementary Figure 11"
	@echo "  make figs12       - Generate Supplementary Figure 12"
	@echo "  make figs13       - Generate Supplementary Figure 13"
	@echo "  make clean        - Remove all generated figures"
	@echo "  make clean-plots  - Remove plots/ directory"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make              # Generate all figures"
	@echo "  make fig1 fig2    # Generate Figure 1 and 2"
	@echo "  make main         # Generate only main figures"
	@echo "  make clean        # Remove generated figures"
