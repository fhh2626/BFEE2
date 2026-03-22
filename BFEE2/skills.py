import json
import re


class AISkillExecutor:
    """Execute AI-generated skill JSON with a strict whitelist."""

    SCHEMA_NAME = "skills-v1"

    def __init__(self, dialog):
        self.dialog = dialog

    @property
    def parent(self):
        return self.dialog.parent

    def _append_exec_message(self, message: str):
        self.dialog._appendExecMessage(message)

    def _coerce_bool(self, value, field_name: str):
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            normalized = value.strip().lower()
            if normalized in ["true", "yes", "1", "on"]:
                return True
            if normalized in ["false", "no", "0", "off"]:
                return False
        raise ValueError(f"{field_name} must be a boolean value")

    def _extract_last_json_code_block(self, ai_text: str):
        matches = re.findall(
            r"```json\s*(\{.*?\})\s*```", ai_text, flags=re.IGNORECASE | re.DOTALL
        )
        if not matches:
            return None
        return matches[-1]

    def _parse_skill_payload(self, ai_text: str):
        json_block = self._extract_last_json_code_block(ai_text)
        if json_block is None:
            return None

        try:
            payload = json.loads(json_block)
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid JSON code block: {exc.msg}") from exc

        if not isinstance(payload, dict):
            raise ValueError("skill payload must be a JSON object")

        if payload.get("bfee_schema") != self.SCHEMA_NAME:
            raise ValueError(
                f'invalid bfee_schema: expected "{self.SCHEMA_NAME}"'
            )

        unexpected_keys = set(payload.keys()) - {"bfee_schema", "skills"}
        if unexpected_keys:
            unexpected = ", ".join(sorted(unexpected_keys))
            raise ValueError(f"unsupported payload keys: {unexpected}")

        return payload

    def _set_ligand_flexibility(self, ligand_type: str):
        ligand_type = ligand_type.strip().lower()
        if ligand_type not in ["flexible", "rigid"]:
            raise ValueError('ligand_type must be "flexible" or "rigid"')

        enabled = ligand_type == "flexible"
        self.parent.geometricAdvancedSettings.considerRMSDCVCheckbox.setChecked(enabled)
        self.parent.alchemicalAdvancedSettings.considerRMSDCVCheckbox.setChecked(
            enabled
        )

    def _skill_protein_protein_geometric(self, args: dict):
        self.parent.selectMDEngineCombobox.setCurrentText("NAMD")
        self.parent.selectStrategyCombobox.setCurrentText("Geometrical")
        self.parent.geometricAdvancedSettings.stratificationRLineEdit.setText("5")
        self.parent.geometricAdvancedSettings.useCUDASOAIntegrator.setChecked(False)
        self.parent.geometricAdvancedSettings.considerRMSDCVCheckbox.setChecked(False)
        self.parent.alchemicalAdvancedSettings.considerRMSDCVCheckbox.setChecked(False)
        self.parent.geometricAdvancedSettings.useGaWTMCheckbox.setChecked(True)

    def _skill_protein_ligand_geometric(self, args: dict):
        ligand_type = str(args.get("ligand_type", "flexible"))
        self._set_ligand_flexibility(ligand_type)

        self.parent.selectMDEngineCombobox.setCurrentText("NAMD")
        self.parent.selectStrategyCombobox.setCurrentText("Geometrical")
        self.parent.geometricAdvancedSettings.stratificationRMSDBoundLineEdit.setText(
            "3"
        )
        self.parent.geometricAdvancedSettings.stratificationRMSDUnboundLineEdit.setText(
            "3"
        )
        self.parent.geometricAdvancedSettings.stratificationRLineEdit.setText("5")
        self.parent.geometricAdvancedSettings.useCUDASOAIntegrator.setChecked(True)
        self.parent.geometricAdvancedSettings.useGaWTMCheckbox.setChecked(False)

    def _skill_protein_ligand_alchemical(self, args: dict):
        ligand_type = str(args.get("ligand_type", "flexible"))
        self._set_ligand_flexibility(ligand_type)

        self.parent.selectMDEngineCombobox.setCurrentText("NAMD")
        self.parent.selectStrategyCombobox.setCurrentText("Alchemical")
        self.parent.alchemicalAdvancedSettings.boundLigandLineEdit.setText("200")
        self.parent.alchemicalAdvancedSettings.boundRestraintsLineEdit.setText("200")
        self.parent.alchemicalAdvancedSettings.unboundLigandLineEdit.setText("100")
        self.parent.alchemicalAdvancedSettings.unboundRestraintsLineEdit.setText("100")
        self.parent.alchemicalAdvancedSettings.doubleWideCheckbox.setChecked(True)
        self.parent.alchemicalAdvancedSettings.useCUDASOAIntegrator.setChecked(True)
        self.parent.alchemicalAdvancedSettings.reEqCheckbox.setChecked(True)
        self.parent.alchemicalAdvancedSettings.LDDMCheckbox.setChecked(False)
        self.parent.alchemicalAdvancedSettings.useWTMLambdaABFCheckbox.setChecked(False)

    def _skill_protein_ligand_lddm(self, args: dict):
        self.parent.selectMDEngineCombobox.setCurrentText("NAMD")
        self.parent.selectStrategyCombobox.setCurrentText("Alchemical")
        self.parent.alchemicalAdvancedSettings.pinDownProCheckbox.setChecked(True)
        self.parent.alchemicalAdvancedSettings.useCUDASOAIntegrator.setChecked(True)
        self.parent.alchemicalAdvancedSettings.LDDMCheckbox.setChecked(True)
        self.parent.alchemicalAdvancedSettings.boundLigandLineEdit.setText("200")
        self.parent.alchemicalAdvancedSettings.unboundLigandLineEdit.setText("100")
        self.parent.alchemicalAdvancedSettings.timestepLineEdit.setText("2.0")
        self.parent.geometricAdvancedSettings.considerRMSDCVCheckbox.setChecked(False)
        self.parent.alchemicalAdvancedSettings.considerRMSDCVCheckbox.setChecked(False)

    def _skill_set_common_fields(self, args: dict):
        if "select_protein" in args:
            self.parent.selectProteinLineEdit.setText(str(args["select_protein"]))
        if "select_ligand" in args:
            self.parent.selectLigandLineEdit.setText(str(args["select_ligand"]))
        if "temperature" in args:
            self.parent.temperatureLineEdit.setText(str(args["temperature"]))

    def _parse_override_args(self, args: dict):
        normalized = {}
        allowed_choices = ["NaCl", "KCl", "CaCl2", "None"]

        for key, value in args.items():
            if key == "consider_rmsd_cv":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "use_gawtm":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "use_cudasoa":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "pin_down_protein":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "use_quaternion_cv":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "reflecting_boundary":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "double_wide":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "re_equilibration":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "use_wtm_lambda_abf":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "use_lddm":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "membrane_protein":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "opls_mixing_rule":
                normalized[key] = self._coerce_bool(value, key)
            elif key == "timestep":
                normalized[key] = str(value)
            elif key == "parallel_runs":
                normalized[key] = str(value)
            elif key == "neutralize_ligand_only":
                normalized[key] = str(value)
                if normalized[key] not in allowed_choices:
                    raise ValueError(
                        "neutralize_ligand_only must be one of NaCl, KCl, CaCl2, None"
                    )
            elif key == "force_field":
                force_field = str(value).strip().lower()
                allowed_force_fields = {
                    "charmm": "CHARMM",
                    "amber": "Amber",
                }
                if force_field not in allowed_force_fields:
                    raise ValueError("force_field must be either CHARMM or Amber")
                normalized[key] = allowed_force_fields[force_field]
            else:
                raise ValueError(f"Unsupported override key: {key}")

        conflicts = []
        if normalized.get("use_gawtm") is True and normalized.get("use_cudasoa") is True:
            conflicts.append("use_gawtm=true conflicts with use_cudasoa=true")
        if normalized.get("use_lddm") is True:
            if normalized.get("consider_rmsd_cv") is True:
                conflicts.append("use_lddm=true conflicts with consider_rmsd_cv=true")
            if normalized.get("use_wtm_lambda_abf") is True:
                conflicts.append(
                    "use_lddm=true conflicts with use_wtm_lambda_abf=true"
                )
            if normalized.get("double_wide") is False:
                conflicts.append("use_lddm=true conflicts with double_wide=false")
            if normalized.get("re_equilibration") is False:
                conflicts.append(
                    "use_lddm=true conflicts with re_equilibration=false"
                )

        if conflicts:
            raise ValueError("; ".join(conflicts))

        return normalized

    def _skill_apply_overrides(self, args: dict):
        normalized = self._parse_override_args(args)

        handlers = {
            "consider_rmsd_cv": lambda v: [
                self.parent.geometricAdvancedSettings.considerRMSDCVCheckbox.setChecked(v),
                self.parent.alchemicalAdvancedSettings.considerRMSDCVCheckbox.setChecked(v),
            ],
            "use_gawtm": lambda v: (
                self.parent.geometricAdvancedSettings.useGaWTMCheckbox.setChecked(v)
            ),
            "use_cudasoa": lambda v: [
                self.parent.geometricAdvancedSettings.useCUDASOAIntegrator.setChecked(v),
                self.parent.alchemicalAdvancedSettings.useCUDASOAIntegrator.setChecked(v),
            ],
            "pin_down_protein": lambda v: [
                self.parent.geometricAdvancedSettings.pinDownProCheckbox.setChecked(v),
                self.parent.alchemicalAdvancedSettings.pinDownProCheckbox.setChecked(v),
            ],
            "use_quaternion_cv": lambda v: [
                self.parent.geometricAdvancedSettings.useOldCvCheckbox.setChecked(v),
                self.parent.alchemicalAdvancedSettings.useOldCvCheckbox.setChecked(v),
            ],
            "reflecting_boundary": lambda v: (
                self.parent.geometricAdvancedSettings.reflectingBoundaryCheckbox.setChecked(v)
            ),
            "double_wide": lambda v: (
                self.parent.alchemicalAdvancedSettings.doubleWideCheckbox.setChecked(v)
            ),
            "re_equilibration": lambda v: (
                self.parent.alchemicalAdvancedSettings.reEqCheckbox.setChecked(v)
            ),
            "use_wtm_lambda_abf": lambda v: (
                self.parent.alchemicalAdvancedSettings.useWTMLambdaABFCheckbox.setChecked(v)
            ),
            "use_lddm": lambda v: (
                self.parent.alchemicalAdvancedSettings.LDDMCheckbox.setChecked(v)
            ),
            "membrane_protein": lambda v: [
                self.parent.geometricAdvancedSettings.memProCheckbox.setChecked(v),
                self.parent.alchemicalAdvancedSettings.memProCheckbox.setChecked(v),
            ],
            "opls_mixing_rule": lambda v: [
                self.parent.geometricAdvancedSettings.OPLSMixingRuleCheckbox.setChecked(v),
                self.parent.alchemicalAdvancedSettings.OPLSMixingRuleCheckbox.setChecked(v),
            ],
            "timestep": lambda v: [
                self.parent.geometricAdvancedSettings.timestepLineEdit.setText(v),
                self.parent.alchemicalAdvancedSettings.timestepLineEdit.setText(v),
            ],
            "parallel_runs": lambda v: (
                self.parent.geometricAdvancedSettings.parallelRunsLineEdit.setText(v)
            ),
            "neutralize_ligand_only": lambda v: [
                self.parent.geometricAdvancedSettings.neutralizeLigOnlyCombobox.setCurrentText(
                    v
                ),
                self.parent.alchemicalAdvancedSettings.neutralizeLigOnlyCombobox.setCurrentText(
                    v
                ),
            ],
            "force_field": lambda v: [
                self.parent.forceFieldCombobox.setCurrentText(v),
                self.parent._changeFFButtonState(),
            ],
        }

        for key, value in normalized.items():
            handlers[key](value)

    def _get_skill_dispatch(self):
        return {
            "protein_protein_geometric": self._skill_protein_protein_geometric,
            "protein_ligand_geometric": self._skill_protein_ligand_geometric,
            "protein_ligand_alchemical": self._skill_protein_ligand_alchemical,
            "protein_ligand_lddm": self._skill_protein_ligand_lddm,
            "set_common_fields": self._skill_set_common_fields,
            "apply_overrides": self._skill_apply_overrides,
        }

    def execute(self, ai_text: str):
        try:
            payload = self._parse_skill_payload(ai_text)
        except ValueError as exc:
            self._append_exec_message(f"Exec: invalid skills payload -> {exc}")
            return

        if payload is None:
            return

        skills = payload.get("skills")
        if not isinstance(skills, list):
            self._append_exec_message(
                'Exec: invalid skills payload -> expected a "skills" array'
            )
            return

        dispatch = self._get_skill_dispatch()
        for item in skills:
            if not isinstance(item, dict):
                self._append_exec_message(
                    "Exec: skip invalid skill entry (must be an object)"
                )
                continue

            name = item.get("name") or item.get("skill")
            args = item.get("args", {})

            if not isinstance(name, str) or name == "":
                self._append_exec_message(
                    "Exec: skip invalid skill entry (missing name)"
                )
                continue
            if not isinstance(args, dict):
                self._append_exec_message(
                    f"Exec: {name} failed: args must be a JSON object"
                )
                continue
            if name not in dispatch:
                self._append_exec_message(f"Exec: skip unsupported skill: {name}")
                continue

            try:
                dispatch[name](args)
                self._append_exec_message(f"Exec: skill {name} -> OK")
            except Exception as e:
                self._append_exec_message(f"Exec: skill {name} failed: {e}")
