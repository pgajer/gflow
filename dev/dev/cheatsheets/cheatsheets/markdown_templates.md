# Markdown Templates

A collection of reusable Markdown structures (templates) for documentation, READMEs, notes, and CRAN submission materials.

---

## 1. Headings

```markdown
# H1: Top-Level Title
## H2: Section Title
### H3: Subsection
#### H4: Sub-subsection
```

---

## 2. Lists

### Unordered list
```markdown
- Item 1
- Item 2
  - Sub-item A
  - Sub-item B
```

### Ordered list
```markdown
1. Step one
2. Step two
   1. Sub-step
   2. Sub-step
3. Step three
```

### Task list
```markdown
- [ ] To do item
- [x] Completed item
```

---

## 3. Code Blocks

### Inline
```markdown
Use `devtools::load_all()` to load the package.
```

### Fenced block
````
```r
devtools::load_all()
Rcpp::compileAttributes()
```
````

Language options: `r`, `bash`, `cpp`, `python`, `make`, etc.

---

## 4. Blockquotes & Callouts

```markdown
> This is a blockquote.
> It can span multiple lines.
```

Custom callout style (common in docs):
```markdown
> **Note:** Be sure to run `Rcpp::compileAttributes()` before building.
> 
> **Warning:** Do not commit compiled `.o` or `.so` files!
```

---

## 5. Tables

```markdown
| Column 1 | Column 2 | Column 3 |
|----------|----------|----------|
| A        | B        | C        |
| 1        | 2        | 3        |
```

---

## 6. Links & Images

### Links
```markdown
[Link text](https://cran.r-project.org)
```

### Images
```markdown
![Alt text](path/to/image.png)
```

---

## 7. Horizontal Rule

```markdown
---
```

---

## 8. Checklists for CRAN / To-do

```markdown
## CRAN Submission Checklist

- [ ] Run `R CMD check --as-cran`
- [ ] Run `rhub::check_for_cran()`
- [ ] Update `cran_comments.md`
- [ ] Confirm DESCRIPTION/NEWS/README are current
- [ ] Verify no `TODO` markers left in code
```

---

## 9. Collapsible Sections (GitHub-flavored)

```markdown
<details>
  <summary>Click to expand</summary>

Content goes here (can include lists, code, etc.)

</details>
```

---

## 10. Badges (GitHub/CRAN READMEs)

```markdown
[![CRAN Status](https://www.r-pkg.org/badges/version/gflow)](https://cran.r-project.org/package=gflow)
[![R-CMD-check](https://github.com/USER/REPO/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/USER/REPO/actions)
```

---

## 11. Footnotes

```markdown
This sentence has a footnote.[^1]

[^1]: Footnote text here.
```

---

## 12. Math (when supported, e.g. GitHub + MathJax)

Inline:
```markdown
Euler's identity: $e^{i\pi} + 1 = 0$
```

Block:
```markdown
$$
\nabla \cdot \vec{E} = \frac{\rho}{\varepsilon_0}
$$
```

---

## 13. Callouts for R Notebooks (R Markdown style)

```markdown
::: note
This is a note callout.
:::

::: warning
This is a warning callout.
:::
```

---

## 14. File Trees

```markdown
project/
├── src/
│   ├── file1.cpp
│   └── file2.cpp
├── R/
└── README.md
```

---

## 15. Example Template for Function Documentation

```markdown
## Function Name

**Description:** Brief description.

**Usage:**
```r
function_name(arg1, arg2, ...)
```

**Arguments:**
- `arg1`: description
- `arg2`: description

**Returns:** What the function outputs.

**Examples:**
```r
function_name(1, 2)
```
```

---
