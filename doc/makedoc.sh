python3 makeargparse.py ../phydget.py prog-desc.md
cat intro.md prog-desc.md footer.md > phydget.md
pandoc -f markdown -t html -o ../phydget.html --standalone --toc --toc-depth=2 -c mdmanstyle1.css --self-contained phydget.md
#pandoc -f markdown -t latex -o phydget.tex --standalone --toc --toc-depth=2 -c mdmanstyle1.css --self-contained phydget.md
pandoc -f markdown -t latex -o ../phydget.pdf --standalone --toc --toc-depth=2 -c mdmanstyle1.css --self-contained phydget.md

